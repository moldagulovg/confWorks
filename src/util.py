from rdkit import Chem
from rdkit.Chem import AllChem
import os 
import time
import subprocess
import shutil
import tempfile
import numpy as np
from tqdm import tqdm

def empty_directory(directory_path):
    try:
        files = os.listdir(directory_path)
        for file in files:
            file_path = os.path.join(directory_path, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
        print("All files deleted successfully.")
    except OSError:
        print("Error occurred while deleting files.")

def copy_conformer_coordinates(source_mol, target_mol, source_confId=-1, target_confId=-1):
    # Check that both molecules have the same number of atoms
    if source_mol.GetNumAtoms() != target_mol.GetNumAtoms():
        raise ValueError("The source and target molecules must have the same number of atoms.")
    
    # Extract the conformer from the source molecule
    source_conf = source_mol.GetConformer(source_confId)
    
    # Create a new conformer for the target molecule (or replace existing one)
    target_conf = target_mol.GetConformer(target_confId)
    
    # Copy coordinates from source conformer to target conformer
    for i in range(source_mol.GetNumAtoms()):
        pos = source_conf.GetAtomPosition(i)
        target_conf.SetAtomPosition(i, pos)
        
    return target_mol

def extract_energy_from_sdf(sdf_file):
    try:
        with open(sdf_file, "r") as file:
            for line in file:
                if line.startswith(" energy:"):
                    # Extract the energy value from the line
                    energy = line.split()[1] # Second element is the energy
                    energy = float(energy.strip())
                    # print(f"Energy: {energy} Hartree")
                    return energy
        # print("Energy not found in the file.")
        return None
    except Exception as e:
        # print(f"An error occurred: {e}")
        return None

def generate_constraints(freeze_atoms, filename="constraints.inp"):
    """
    Generates an XTB constraints.inp file to freeze specified atoms.

    Parameters:
    - freeze_atoms (list of int): List of atom indices (0-based) to freeze.
    - filename (str): Output constraint filename (default: constraints.inp).
    """

    # Convert to 1-based indexing for xtb
    freeze_atoms_xtb = [i + 1 for i in freeze_atoms]

    # Write constraint file
    with open(filename, "w") as f:
        f.write("$fix\n")
        f.write("atoms: " + ", ".join(map(str, freeze_atoms_xtb)) + "\n")
        f.write("$end\n")

    # print(f"Constraints file '{filename}' created successfully."



# TODO make this function work silently, without much verbosity
# TODO make RMSD filter, if multiple conformers are close enough delete one of them and keep the other.
def xtb_SP(mol_object,
            gfn_xtb=2,
            solvent=None,
            unpaired_e=None,
            verbose=False,
            silent=False,
            confId=False,
            charge=0,
            freeze_atoms=None
            ):
    
    #############################################################################
    #############################################################################

    xtb_flags = []

    if gfn_xtb==2 or gfn_xtb==1 or gfn_xtb==0:
        xtb_flags+=["--gfn", str(gfn_xtb)]
    elif gfn_xtb=='ff':
        xtb_flags+=["--gfnff",]
    else:
        print("Invalid xTB parametrization. Choose from 2 (default), 1 or 0. Integer input.")
    
    if isinstance(charge, int):
        xtb_flags += ["--chrg", str(charge)]
    else:
        raise ValueError(f'charge shall be integer, {type(charge)} was given instead')

    if solvent is not None:
        ##check if the solvent is correct
        available_solvents = ["acetone", 'acetonitrile', 'benzene', 'ch2cl2', 'chcl3', 'cs2', 
                                    'dioxane', 'dmf', 'dmso', 'ether', 'ethylacetate', 'furane', 
                                    'hexane', 'methanol', 'nitromethane', 'toluene', 'thf', 'water']
        if not isinstance(solvent, str):
            raise ValueError(f"Solvent input shall be in string format, but {type(solvent)} was given instead.")
        
        elif solvent.lower() not in available_solvents:
            raise ValueError(f"Invalid solvent choice. \nPlease choose one from the list: ", available_solvents)

        xtb_flags+=["--alpb", solvent]

    if unpaired_e is not None:
        ##check if the format is correct
        xtb_flags+=["--uhf", str(unpaired_e)]
    
    if verbose is True:
        xtb_flags+=['--verbose',]
    elif silent is True:
        xtb_flags+=['--silent',]
    else:
        pass

    

    path = os.getcwd()
    numConfs = mol_object.GetNumConformers()

    # TODO add confId value checker, should be either False or int value within numConfs.

    for confX in tqdm(range(numConfs)):
        if confId:
            if confX != confId:
                continue
        try:
            temp_dir = tempfile.mkdtemp(prefix='temp_geom_opt_')
            os.chdir(temp_dir)
            
            Chem.SDWriter('out.sdf').write(mol_object, confId=confX)

            if freeze_atoms is not None:
                generate_constraints(freeze_atoms, filename="constraints.inp")
                xtb_flags+=['--input', 'constraints.inp']
            
            output_filename = "xtb_output.txt"

            start = time.time()
            command = ["xtb", "out.sdf"] + xtb_flags + [">>", f"{output_filename}", ]
            result = subprocess.run(command, 
                                    capture_output=True,
                                    text=True,
                                    check=True  # Raises CalledProcessError if the command returns a non-zero exit code
                                    )

            # The captured standard output is in 'result.stdout'
            output_string = result.stdout

            # Process each line of the output string
            for line in output_string.splitlines():
                if "TOTAL ENERGY" in line:
                    words = line.split()
                    # The energy value is the second to last element
                    total_energy = float(words[-3])
                
            # print(temp_dir)
            # print(" ".join(command))
            end = time.time()  
            # print(f'elapsed time for optimization: {end-start}')

            # except FileNotFoundError:
            #     print(f"Error: The file '{output_filename}' was not found.")
            #     return None
            # except (ValueError, IndexError):
            #     print("Error: Could not parse the total energy from the file.")
            #     return None
            energy = total_energy
            conf = mol_object.GetConformer(confX)
            conf.SetIntProp("conf_id", confX)
            conf.SetDoubleProp("conf_energy", energy)

        finally:
            os.chdir(path)
            shutil.rmtree(temp_dir)

    return mol_object


def optimize_molecule(mol_object,
                      gfn_xtb=2,
                      opt_lvl='normal',
                      num_opt_cycles=None,
                      solvent=None,
                      unpaired_e=None,
                      verbose=False,
                      silent=False,
                      confId=False,
                      charge=0,
                      freeze_atoms=None
                      ):
    
    #############################################################################
    #############################################################################

    xtb_flags = []

    if gfn_xtb==2 or gfn_xtb==1 or gfn_xtb==0:
        xtb_flags+=["--gfn", str(gfn_xtb)]
    elif gfn_xtb=='ff':
        xtb_flags+=["--gfnff",]
    else:
        print("Invalid xTB parametrization. Choose from 2 (default), 1 or 0. Integer input.")
    
    if isinstance(charge, int):
        xtb_flags += ["--chrg", str(charge)]
    else:
        raise ValueError(f'charge shall be integer, {type(charge)} was given instead')

    # check if chosen opt lvl is valid
    valid_opt_lvls=["crude", "sloppy", "loose", "lax", "normal", "tight", "vtight", "extreme"]
    if opt_lvl not in valid_opt_lvls:
        raise ValueError(f'Chosen xTB geometry optimization level is not valid.\n'+
              'Please choose one from the list: ', valid_opt_lvls)

    xtb_flags+=["--opt", opt_lvl]

    if num_opt_cycles is not None:
        if not isinstance(num_opt_cycles, int):
            raise ValueError(f"Num_opt_cycles shall be in int format, but {type(num_opt_cycles)} was given instead.")

        xtb_flags+=["--cycles", str(num_opt_cycles)]

    if solvent is not None:
        ##check if the solvent is correct
        available_solvents = ["acetone", 'acetonitrile', 'benzene', 'ch2cl2', 'chcl3', 'cs2', 
                                    'dioxane', 'dmf', 'dmso', 'ether', 'ethylacetate', 'furane', 
                                    'hexane', 'methanol', 'nitromethane', 'toluene', 'thf', 'water']
        if not isinstance(solvent, str):
            raise ValueError(f"Solvent input shall be in string format, but {type(solvent)} was given instead.")
        
        elif solvent.lower() not in available_solvents:
            raise ValueError(f"Invalid solvent choice. \nPlease choose one from the list: ", available_solvents)

        xtb_flags+=["--alpb", solvent]

    if unpaired_e is not None:
        ##check if the format is correct
        xtb_flags+=["--uhf", str(unpaired_e)]
    
    if verbose is True:
        xtb_flags+=['--verbose',]
    elif silent is True:
        xtb_flags+=['--silent',]
    else:
        pass

    

    path = os.getcwd()
    numConfs = mol_object.GetNumConformers()

    # TODO add confId value checker, should be either False or int value within numConfs.

    for confX in range(numConfs):
        if confId:
            if confX != confId:
                continue
        try:
            temp_dir = tempfile.mkdtemp(prefix='temp_geom_opt_')
            os.chdir(temp_dir)
            
            Chem.SDWriter('out.sdf').write(mol_object, confId=confX)

            if freeze_atoms is not None:
                generate_constraints(freeze_atoms, filename="constraints.inp")
                xtb_flags+=['--input', 'constraints.inp']
            start = time.time()
            subprocess.run(["xtb", "out.sdf"] + xtb_flags)
            end = time.time()  
            print(f'elapsed time for optimization: {end-start}')

            # TODO add checking if xtb optimization crashed for the file or was successfull
            # TODO get geometry optimization statistics including (1) elapsed time for each input, (2) num of iterations for each, (3) num successful, (4) num failed.
            energy = extract_energy_from_sdf('xtbopt.sdf')
            suppl = Chem.ForwardSDMolSupplier('xtbopt.sdf', sanitize=False)
            for source_mol in suppl:
                # Example usage:
                # Assuming source_mol and target_mol are already defined and have conformers
                mol_object = copy_conformer_coordinates(source_mol, mol_object, source_confId=-1, target_confId=confX)
            conf = mol_object.GetConformer(confX)
            conf.SetIntProp("conf_id", confX)
            conf.SetDoubleProp("conf_energy", energy)

        finally:
            os.chdir(path)
            shutil.rmtree(temp_dir)  

    return mol_object


def conformer_search(mol, 
                    gfn_xtb=None, 
                    output_dir=None, 
                    solvent=None, 
                    threads=None, 
                    topo_change=False,
                    settings='normal',
                    charge=0,
                    ):
    
    crest_flags = []
    
    if isinstance(charge, int):
        crest_flags += ["--chrg", str(charge)]
    else:
        raise ValueError(f'charge shall be integer, {type(charge)} was given instead')
    
    if gfn_xtb==2 or gfn_xtb==1 or gfn_xtb==0:
        crest_flags+=["--gfn",str(gfn_xtb)]
    elif gfn_xtb=='ff':
        crest_flags+=["--gfnff",]
    else:
        raise ValueError("Invalid xTB parametrization. Choose from 2 (default), 1 or 0. Integer input.")
    if threads != None:
        crest_flags+=['-T',str(threads)]
        
    if solvent is not None:
        ##check if the solvent is correct
        available_solvents = ["acetone", 'acetonitrile', 'benzene', 'ch2cl2', 'chcl3', 'cs2', 
                                    'dioxane', 'dmf', 'dmso', 'ether', 'ethylacetate', 'furane', 
                                    'hexane', 'methanol', 'nitromethane', 'toluene', 'thf', 'water']
        if not isinstance(solvent, str):
            raise ValueError(f"Solvent input shall be in string format, but {type(solvent)} was given instead.")
        elif solvent.lower() not in available_solvents:
            raise ValueError(f"Invalid solvent choice. \nPlease choose one from the list: ", available_solvents)
        crest_flags+=["--alpb", solvent]
        
    if topo_change == True:
            crest_flags+=['--noreftopo']
    else:
        pass
    
    settings_list = ['quick', 'squick', 'mquick']
    if settings in settings_list:
        crest_flags+=[f'--{settings}']
    elif settings == 'normal':
        pass
    else:
        raise ValueError(f"Invalid speed setting selected. \nPlease choose one from the list: ", settings_list)
        

    cwd_path = os.getcwd()
    try:
        temp_dir = tempfile.mkdtemp(prefix='temp_geom_opt_')
        os.chdir(temp_dir)
        writer = Chem.SDWriter("input.sdf")
        writer.write(mol)
        writer.close()
    
        start = time.time()
        process = subprocess.run(['crest', f'{"input.sdf"}'] + crest_flags,)
        end = time.time()  
        print(f'elapsed time for conformer search: {end-start}')

        # Parse the output SDF from CREST
        output_sdf = os.path.join(temp_dir, "crest_conformers.sdf")
        if not os.path.exists(output_sdf):
            raise FileNotFoundError(f"CREST output file not found: {output_sdf}")
        
        # Load conformers from SDF into RDKit mol object
        conformer_supplier = Chem.SDMolSupplier(output_sdf, removeHs=False)
        base_mol = Chem.Mol(conformer_supplier[0])
        
        # TODO consider loading energies of each conformers
        for i, mol in enumerate(conformer_supplier):
            if i != 0:
                conf = mol.GetConformer()
                base_mol.AddConformer(conf)

        print(f'Successfully sampled {base_mol.GetNumConformers()} conformers')
        # Extract energies from the second line of each conformer in the SDF file
        conformer_energies = {}
        with open(output_sdf, 'r') as sdf_file:
            lines = sdf_file.readlines()
        
        i = 0  # index of the conformer
        for line_index, line in enumerate(lines):
            # TODO below line might be incorrect, probably it should be 'energy: '
            if "Energy =" in line:  # Check for the line with energy info
                energy_str = line.split("=")[-1].strip()  # Extract the energy part
                energy = float(energy_str)
                conformer_energies[i] = energy
                i += 1  # Move to the next conformer
    
    finally:
        os.chdir(cwd_path)
        if output_dir:
            # Ensure output directory doesn't exist
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            shutil.copytree(temp_dir, output_dir)
        shutil.rmtree(temp_dir)

    return base_mol, conformer_energies


def get_elements_coordinates(mol, confId=-1):
    pt = Chem.GetPeriodicTable()
    # Get the atomic numbers and map them to element symbols
    elements = np.array([pt.GetElementSymbol(atom.GetAtomicNum()) for atom in mol.GetAtoms()])

    # Get the conformer with 3D coordinates
    conformer = mol.GetConformer(id=confId)

    # Extract the 3D coordinates for each atom
    coordinates = np.array([conformer.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    coordinates = coordinates.astype(float)  # Convert to numpy float array
    return elements, coordinates


def calc_rmsd(mol, correct_mol):
    '''
    https://github.com/charnley/rmsd.git

    Function takes two RDKit Mol objects with embedded conformers.
    
    For each conformer of mol_obj calculates RMSD against all conformers of correct_mol.

    Returns list of lists of size N_conf_mol_obj x N_conf_correct_mol.
    '''
    correct_mol = Chem.RemoveAllHs(correct_mol, sanitize=False)
    
    path = os.getcwd()
    rmsd_values_for_conformer = []
    try:
        temp_dir = tempfile.mkdtemp(prefix='temp_geom_opt_')
        os.chdir(temp_dir)

        values = []

        # use temporary files directory for working with XYZ files
        mol = Chem.RemoveAllHs(mol, sanitize=False)
        for i in range(mol.GetNumConformers()):
            Chem.rdmolfiles.MolToXYZFile(mol, "mol.xyz", confId=i)
            for j in range(correct_mol.GetNumConformers()):
                Chem.rdmolfiles.MolToXYZFile(correct_mol, "correct.xyz", confId=j)
                command = f"calculate_rmsd  --reorder --no-hydrogen correct.xyz mol.xyz"

                # Run the command and capture the output
                value = subprocess.check_output(command, shell=True, text=True)
                values.append(value)
        rmsd_values_for_conformer.append([float(value.replace('\n', '')) for value in values])
    finally:
        os.chdir(path)
        shutil.rmtree(temp_dir)  
    return rmsd_values_for_conformer

def calc_rmsd_single_mol(mol, verbose = False):
    '''
    https://github.com/charnley/rmsd.git

    Function takes an RDKit Mol object with embedded conformers.
    
    For each conformer of mol_obj calculates RMSD against all other conformers.

    Returns list of lists of size (N_conf x (N_conf - 1)) / 2.
    '''
    mol = Chem.RemoveAllHs(mol, sanitize=False)
    
    path = os.getcwd()
    rmsd_values_for_conformer = []
    try:
        temp_dir = tempfile.mkdtemp(prefix='temp_geom_opt_')
        os.chdir(temp_dir)

        values = []

        # use temporary files directory for working with XYZ files
        mol = Chem.RemoveAllHs(mol, sanitize=False)
        for i in tqdm(range(mol.GetNumConformers())):
            Chem.rdmolfiles.MolToXYZFile(mol, "mol.xyz", confId=i)
            for j in range(i):
                Chem.rdmolfiles.MolToXYZFile(mol, "correct.xyz", confId=j)
                command = f"calculate_rmsd  --reorder --no-hydrogen correct.xyz mol.xyz"

                # Run the command and capture the output
                value = subprocess.check_output(command, shell=True, text=True)
                values.append(value)
        rmsd_values_for_conformer.append([float(value.replace('\n', '')) for value in values])
    finally:
        os.chdir(path)
        shutil.rmtree(temp_dir)  
    return rmsd_values_for_conformer

def calc_rmsd_matrix_single_mol(mol, verbose = False):
    '''
    https://github.com/charnley/rmsd.git

    Function takes an RDKit Mol object with embedded conformers.
    
    For each conformer of mol_obj calculates RMSD against all other conformers.

    Returns list of lists of size (N_conf x (N_conf - 1)) / 2.
    '''
    mol = Chem.RemoveAllHs(mol, sanitize=False)
    numConfs = mol.GetNumConformers()
    path = os.getcwd()
    rmsd_matrix = np.zeros((numConfs, numConfs))
    try:
        temp_dir = tempfile.mkdtemp(prefix='temp_geom_opt_')
        os.chdir(temp_dir)

        # use temporary files directory for working with XYZ files
        mol = Chem.RemoveAllHs(mol, sanitize=False)
        for i in tqdm(range(numConfs)):
            Chem.rdmolfiles.MolToXYZFile(mol, "mol.xyz", confId=i)
            for j in range(i):
                Chem.rdmolfiles.MolToXYZFile(mol, "correct.xyz", confId=j)
                command = f"calculate_rmsd  --reorder --no-hydrogen correct.xyz mol.xyz"

                # Run the command and capture the output
                value = subprocess.check_output(command, shell=True, text=True)
                value = float(value.replace('\n', ''))
                rmsd_matrix[i][j], rmsd_matrix[j][i] = value, value
    finally:
        os.chdir(path)
        shutil.rmtree(temp_dir)  
    return rmsd_matrix


def read_multiconf_sdf(filepath):
    suppl = Chem.ForwardSDMolSupplier(filepath, sanitize=False)

    numConfs = 0
    for confX, source_mol in enumerate(suppl):
        numConfs += 1
        if confX == 0:
            mol = source_mol
        else:
            conf = source_mol.GetConformer(-1)
            new_confId = mol.AddConformer(conf, assignId=True)
            if new_confId != confX:
                raise ValueError('different conformer IDs!')
            
        new_conf = mol.GetConformer(confX)
        new_conf.SetDoubleProp('conf_energy', source_mol.GetDoubleProp('conf_energy'))

    if mol.GetNumConformers() == numConfs:
        return mol
    else:
        raise ValueError('Not enough conformers were read.')

def save_multiconf_sdf(mol, filepath):
    mol_copy = Chem.Mol(mol)
    if mol_copy.GetNumConformers() == 0:
        raise ValueError('No conformers in this molecule.')
    
    writer = Chem.SDWriter(filepath)
    for confX in range(mol_copy.GetNumConformers()):
        conf = mol_copy.GetConformer(confX)
        mol_copy.SetIntProp("conf_id", confX)
        if conf.HasProp("conf_energy"):
            energy = conf.GetDoubleProp("conf_energy")
            mol_copy.SetDoubleProp("conf_energy", energy)
        writer.write(mol_copy, confId=confX)
    writer.close()


# TODO calculate SP energies using xtb 

