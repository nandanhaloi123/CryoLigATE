import torch
import numpy as np
import argparse
from pathlib import Path
import sys
import gemmi
import scipy.ndimage
import mrcfile

from CryoLigate.model import SCUNet                

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
TARGET_VOXEL_SIZE = 0.5 
GRID_SIZE = 64

def load_model(weights_path, device):
    """Loads the SCUNet model."""
    print(f"Loading model weights from: {weights_path}")
    model = SCUNet(in_nc=1, window_size=4).to(device)
    state_dict = torch.load(weights_path, map_location=device)
    if list(state_dict.keys())[0].startswith('module.'):
        state_dict = {k[7:]: v for k, v in state_dict.items()}
    model.load_state_dict(state_dict)
    model.eval()
    return model

def save_mrc_with_origin(data, filepath, unit_cell, origin_angstroms):
    """Saves a numpy array as .mrc with correct PHYSICAL origin and DIMENSIONS."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with mrcfile.new(str(filepath), overwrite=True) as mrc:
        mrc.set_data(data.astype(np.float32))
        mrc.header.cella.x = unit_cell.a
        mrc.header.cella.y = unit_cell.b
        mrc.header.cella.z = unit_cell.c
        mrc.voxel_size = unit_cell.a / data.shape[2] 
        mrc.header.origin.x = float(origin_angstroms[0])
        mrc.header.origin.y = float(origin_angstroms[1])
        mrc.header.origin.z = float(origin_angstroms[2])
        mrc.update_header_stats()

def get_cropped_and_resampled_density(gemmi_grid, centroid, grid_size, target_voxel):
    """Extracts the density box around the centroid, mimicking the training prep."""
    box_angstroms = grid_size * target_voxel
    phys_origin = np.array([centroid[0] - box_angstroms/2.0, 
                            centroid[1] - box_angstroms/2.0, 
                            centroid[2] - box_angstroms/2.0], dtype=np.float32)
                            
    x = phys_origin[0] + np.arange(grid_size) * target_voxel
    y = phys_origin[1] + np.arange(grid_size) * target_voxel
    z = phys_origin[2] + np.arange(grid_size) * target_voxel
    xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')
    phys_coords = np.stack([xv.ravel(), yv.ravel(), zv.ravel()], axis=1)
    
    mat = np.array(gemmi_grid.unit_cell.frac.mat.tolist())
    vec = np.array(gemmi_grid.unit_cell.frac.vec.tolist())
    frac_coords = (mat @ phys_coords.T).T + vec
    grid_indices = frac_coords * np.array([gemmi_grid.nu, gemmi_grid.nv, gemmi_grid.nw])
    
    min_idx = np.floor(np.min(grid_indices, axis=0)).astype(int) - 2
    max_idx = np.ceil(np.max(grid_indices, axis=0)).astype(int) + 2
    shape = (max_idx - min_idx + 1).tolist()
    
    sub_arr = np.array(gemmi_grid.get_subarray(min_idx.tolist(), shape))
    local_indices = grid_indices - min_idx
    resampled_1d = scipy.ndimage.map_coordinates(sub_arr, local_indices.T, order=3, mode='grid-wrap')
    
    return resampled_1d.reshape((grid_size, grid_size, grid_size)).astype(np.float32), phys_origin

def run_real_world_inference(model, map_path, pdb_path, resname, resid, output_dir, device):
    print(f"\n--- Processing Complex: {Path(pdb_path).name} ---")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load Structure and Map
    print(f"Reading PDB/CIF: {pdb_path}")
    st = gemmi.read_structure(str(pdb_path))
    
    print(f"Reading Map: {map_path}")
    m = gemmi.read_ccp4_map(str(map_path))
    m.setup(0.0, gemmi.MapSetup.Full)

    # 2. Parse all residues to separate Protein from the specific Target Ligand
    target_residues = []
    protein_residues = []
    protein_atoms = []
    
    # Standardize inputs for comparison
    target_resname = resname.upper()
    target_resid = str(resid)

    for model_obj in st:
        for chain in model_obj:
            for res in chain:
                # Check for the exact Sequence Number AND Residue Name
                if str(res.seqid.num) == target_resid:
                    if res.name.upper() == target_resname:
                        target_residues.append((chain.name, res.seqid.num, res.name, res))
                    else:
                        print(f"WARNING: Found residue {res.name} at Chain {chain.name}:{target_resid}, but expected {target_resname}.")
                elif res.het_flag == 'A': # Amino acids
                    protein_residues.append(res)
                    for atom in res:
                        protein_atoms.append(atom)

    if not target_residues:
        print(f"ERROR: Ligand {target_resname} at Residue {target_resid} not found in the provided structure.")
        return

    # 3. Process the exact ligand instance found
    for (chain_name, res_num, res_name, lig_res) in target_residues:
        identifier = f"{res_name}_Chain{chain_name}_{res_num}"
        print(f"  -> Cropping and inferring for {identifier}")

        # Calculate Centroid
        lig_coords = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_res])
        if len(lig_coords) == 0:
            print(f"     Skipping {identifier} - no atoms found.")
            continue
        centroid = np.mean(lig_coords, axis=0)

        # Extract Cropped Density
        density, phys_origin = get_cropped_and_resampled_density(m.grid, centroid, GRID_SIZE, TARGET_VOXEL_SIZE)
        
        # Save a copy of the raw density before normalization to export later
        raw_density = density.copy()

        # Apply the same normalization used in training for the model
        density = (density - np.mean(density)) / (np.std(density) + 1e-6)

        # Generate Protein Mask (Binary map indicating protein atom presence)
        protein_mask = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=np.float32)
        for atom in protein_atoms:
            pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
            v = np.round((pos - phys_origin) / TARGET_VOXEL_SIZE).astype(np.int32)
            if np.all((v >= 0) & (v < GRID_SIZE)):
                protein_mask[v[0], v[1], v[2]] = 1.0

        # Create Model Input Tensor: Shape (1, 2, D, H, W)
        # input_tensor = np.stack([density, protein_mask], axis=0)
        input_tensor = np.stack([density], axis=0)
        input_torch = torch.from_numpy(input_tensor).float().unsqueeze(0).to(device)

        # Run Inference
        with torch.no_grad():
            pred_tensor = model(input_torch)
            
        # Convert output to saveable format (transpose to Z, Y, X for MRC standard)
        pred_np = pred_tensor.cpu().numpy().squeeze().T

        # --- SAVE OUTPUTS ---
        crop_cell = gemmi.UnitCell(
            GRID_SIZE * TARGET_VOXEL_SIZE,
            GRID_SIZE * TARGET_VOXEL_SIZE,
            GRID_SIZE * TARGET_VOXEL_SIZE,
            90, 90, 90 
        )
        
        # 1. Save Predicted MRC
        out_pred = output_dir / f"refined_{identifier}.mrc"
        save_mrc_with_origin(pred_np, out_pred, crop_cell, phys_origin)
        print(f"     Saved Prediction: {out_pred.name}")

        # 2. Save Original Cropped MRC (Raw, un-normalized)
        out_orig = output_dir / f"original_{identifier}.mrc"
        save_mrc_with_origin(raw_density.T, out_orig, crop_cell, phys_origin)
        print(f"     Saved Original:   {out_orig.name}")

        # 3. Save Local PDB (Ligand + Pocket)
        new_st = gemmi.Structure()
        new_st.cell = st.cell 
        new_model = gemmi.Model("1")
        new_chain = gemmi.Chain("A")
        new_chain.add_residue(lig_res)
        
        pocket_radius_sq = 20.0**2
        for pres in protein_residues:
            if len(pres) > 0:
                # Use first atom to quickly check bounding distance
                atom_pos = pres[0].pos
                dist_sq = (atom_pos.x - centroid[0])**2 + (atom_pos.y - centroid[1])**2 + (atom_pos.z - centroid[2])**2
                if dist_sq < pocket_radius_sq:
                    new_chain.add_residue(pres)
                    
        new_model.add_chain(new_chain)
        new_st.add_model(new_model)
        
        out_pdb = output_dir / f"local_{identifier}.pdb"
        new_st.write_pdb(str(out_pdb))
        print(f"     Saved Local PDB:  {out_pdb.name}")

    print("\nInference Complete!")

def main():
    parser = argparse.ArgumentParser(description="CryoLigate Inference")
    parser.add_argument("--weights", type=str, default="weights/cryoligate_v2.0.0.pth", help="Path to trained weights")
    parser.add_argument("--map", type=str, required=True, help="Full experimental density map (.mrc/.map)")
    parser.add_argument("--pdb", type=str, required=True, help="Full modeled complex (.pdb/.cif)")
    
    # Specific target arguments
    parser.add_argument("--resname", type=str, required=True, help="3-letter code of the ligand (e.g., ATP)")
    parser.add_argument("--resid", type=str, required=True, help="Residue sequence number (e.g., 501)")
    
    parser.add_argument("--outdir", type=str, default=None, help="Directory to save the refined local boxes (default: same folder as --map)")
    
    args = parser.parse_args()
    
    # Default outdir to the map's parent directory if not specified
    if args.outdir is None:
        args.outdir = Path(args.map).parent
    
    for path in [args.weights, args.map, args.pdb]:
        if not Path(path).exists():
            print(f"ERROR: File not found -> {path}")
            sys.exit(1)

    model = load_model(args.weights, DEVICE)
    run_real_world_inference(model, args.map, args.pdb, args.resname, args.resid, args.outdir, DEVICE)

if __name__ == "__main__":
    main()