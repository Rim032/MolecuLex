import time
import os
import re
import csv
import sys

import argparse
from rdkit import Chem
import pubchempy as pcp
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


os.system('')



max_request_size = 256
desired_extension = ".txt"

log_types = [
    "LOG",      #0
    "ERROR",    #1
    "WARNING"   #2
]

display_drug_attributes = True
application_version = 0.97

console_noprint = False
save_csv_data = False
record_full_data = False
disable_stat_report = False

violation_count = 0
likely_drug_count = 0
unlikely_drug_count = 0


csv_file_name = "untitled_molecuLex_data"
csv_data = [
    ["Compound Name", "PubChem ID", "Molecular Weight (g/mol)", "H-Donors", "H-Acceptors", "LogP", "Drug Possibility", "Lipinski Violations"]
]




def print_log(log_type: int, message: str) -> None:
    if message is None or log_type is None:
        return
    print(f"[{log_types[log_type]}] - {message}")


def append_csv_molecule_data(molecule_data: dict) -> None:
    if molecule_data is None:
        print_log(1, f"Cannot append molecule data to CSV data list. Passed data is invalid.")
        return

    global csv_data
    csv_data.append(molecule_data)


def save_csv_molecule_data() -> None:
    try:
        with open(f"{csv_file_name}.csv", mode='w', newline='', encoding='utf-8') as csv_file:
            file_writer = csv.writer(csv_file, delimiter=',')
            file_writer.writerows(csv_data)
        print_log(0, f"Successfully saved data to {csv_file_name}.csv")
    except Exception as err:
        print_log(1, f"Failed to save CSV: {err}")


def record_molecule_info(molecule_data: dict) -> None:
    if console_noprint is True or molecule_data is None:
        return

    mw: float = molecule_data["mw"]
    log_p: float = molecule_data["log_p"]
    h_acceptors: int = molecule_data["h_acceptors"]
    h_donors: int = molecule_data["h_donors"]
    is_drug_likely: bool = molecule_data["is_drug_likely"]
    name: str = molecule_data["name"]
    
    WIDTH: int = 62
    LABEL_W: int = 20
    DATA_W: int = 30
    
    try:
        tag = "[ LIKELY DRUG CANDIDATE ]" if is_drug_likely else "[ UNLIKELY DRUG CANDIDATE ]"
        print(f"\n╔{'═' * (WIDTH-2)}╗")

        molecule_id = molecule_data["id"]
        prefix = f"• CID ({molecule_id}): "
        available_space = (WIDTH - 4) - len(prefix)

        if len(name) > available_space:
            display_name = name[:available_space - 2] + ".."
        else:
            display_name = name
        full_content = f"{prefix}{display_name}"
        
        print(f"║ {full_content:<{WIDTH - 4}} ║")
        print(f"╠{'═' * (WIDTH-2)}╣")

        print(f"║ { 'Molecular Weight:':<{LABEL_W}}{mw:>{DATA_W}.4f} g/mol {' ':<2}║")
        print(f"║ { 'LogP:':<{LABEL_W}}{log_p:>{DATA_W}.4f} {' ':<8}║")
        print(f"║ { 'H-Donors:':<{LABEL_W}}{h_donors:>{DATA_W}} {' ':<8}║")
        print(f"║ { 'H-Acceptors:':<{LABEL_W}}{h_acceptors:>{DATA_W}} {' ':<8}║")

        if display_drug_attributes is True:
            print(f"╟{'─' * (WIDTH-2)}╢")
            print(f"║ {tag:^{WIDTH-4}} ║")
            
            for rule_violator in molecule_data["violations"]:
                clean_violation = (rule_violator[:WIDTH-10] + '..') if len(rule_violator) > (WIDTH-10) else rule_violator
                print(f"║ ! {clean_violation:<{WIDTH-5}}║")
        print(f"╚{'═' * (WIDTH-2)}╝")
    except Exception as err:
        print_log(1, f"Printing issue. {err}.")


def record_molecule_summary(id_list: list[str]) -> None:
    if disable_stat_report is True or id_list is None:
        return
    
    WIDTH: int = 62
    LABEL_W: int = 20
    DATA_W: int = 30
    
    try:
        print("\n\n")
        print(f"\n╔{'═' * (WIDTH-2)}╗")
        full_content = f"+++ SUMMARY REPORT +++"
        
        print(f"║ {full_content:<{WIDTH - 4}} ║")
        print(f"╠{'═' * (WIDTH-2)}╣")

        print(f"║ { 'Total Compounds:':<{LABEL_W}}{len(id_list):>{DATA_W}}\t   {' ':<2}║")
        print(f"║ { 'Total Violations:':<{LABEL_W}}{violation_count:>{DATA_W}}\t   {' ':<2}║")
        print(f"║ { '# Likely Drugs:':<{LABEL_W}}{likely_drug_count:>{DATA_W}}\t   {' ':<2}║")
        print(f"║ { '# Unlikely Drugs:':<{LABEL_W}}{unlikely_drug_count:>{DATA_W}}\t   {' ':<2}║")
        
        print(f"╚{'═' * (WIDTH-2)}╝")
    except Exception as err:
        print_log(1, f"Printing issue. {err}.")


def create_local_id_list(min_id: int, max_id: int) -> list[str]:
    if max_id > sys.maxsize:
        print_log(1, f"Too large of a maximum ID value. MAX INT SIZE: {sys.maxsize}")
        return

    new_id_list: list[str] = []
    for i in range(max_id-min_id):
        new_id_list.append(str(i+min_id+1))

    return new_id_list


def read_id_file(file: str) -> str:
    if file is None or os.path.isfile(file) is False:
        print_log(1, "File is invalid or malformed.")
        return

    file_content: str = ""
    with open(file, "r") as file_reader:
        file_content = file_reader.read()

    return file_content


def format_id_list(ids: str) -> list[str]:
    if ids is None:
        print_log(1, "An invalid list of IDs was passed.")
        return

    id_list: list[str] = re.split(r",|\s+", ids)
    return id_list


def api_bar_print(id_list: list[str], bar_length: int, compound_index: int, compound_id: int) -> None:
    api_progress = compound_index / len(id_list)
    filled_len = int(bar_length * api_progress)
    api_bar = "▓" * filled_len + "░" * (bar_length - filled_len)
    que_percent = int(api_progress * 100)
    sys.stdout.write(f"\rPROGRESS: [{api_bar}] {que_percent}% | Analyzing CID: {compound_id:<10}")
    sys.stdout.flush()


def get_id_smiles(id_list: list[str]) -> list[str]:
    if id_list is None:
        print_log(1, "An invalid list of formatted IDs was passed.")
        return

    smiles_list: list[str] = []
    name_list: list[str] = []

    comp_index: int = 0
    print(f"\n{'—'*20} DATA ACQUISITION PHASE: INITIATED {'—'*20}")
    print_log(0, f"Connecting to PubChem API for {len(id_list)} compounds...")

    try:
        for i in range(0, len(id_list), max_request_size):
            id_list_chunk = id_list[i : i + max_request_size]
            fetched_compounds = pcp.get_compounds(id_list_chunk, 'cid')
            for pub_compound in fetched_compounds:     
                api_bar_print(id_list, 30, comp_index, id_list[comp_index])
                    
                if pub_compound is not None:
                    smiles_list.append(pub_compound.smiles)
                    name_list.append(pub_compound.iupac_name)

                comp_index = comp_index + 1
    except Exception as err:
        print("")
        print_log(1, f"API Issue. {err}.")
    print(f"\n{'—'*20} DATA ACQUISITION PHASE: COMPLETED {'—'*20}")
    
    return smiles_list, name_list


def drug_candidate_check(mw: float, h_donors: int, h_acceptors: int, logp: float) -> bool:
    if mw is None or h_donors is None or h_acceptors is None or logp is None:
        print_log(1, "Cannot check the drug properties of a compound. Invalid arguments were passed.")
        return "N/A Couldn't parse", False
    violations: list[str] = []
    
    if mw > 500:
        violations.append(f"Molecular Weight (g/mol) > 500)")
    if h_donors > 5:
        violations.append(f"Hydrogen bond donors > 5)")
    if h_acceptors > 10:
        violations.append(f"Hydrogen bond acceptors > 10)")
    if logp > 5:
        violations.append(f"Octanol-water partition coefficient (LogP) > 5)")
        
    if len(violations) >= 1:
        global violation_count
        violation_count = violation_count + len(violations)

        global unlikely_drug_count
        unlikely_drug_count = unlikely_drug_count + 1
        
        return violations, False
    global likely_drug_count
    likely_drug_count = likely_drug_count + 1
    return violations, True


def parse_compound_properties(smiles_list: list[str], names_list: list[str], id_list: list[str]) -> None:
    if smiles_list is None:
        print_log(1, "An invalid SMILES list was passed.")
        return
    
    for smiles, name, pc_id in zip(smiles_list, names_list, id_list):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            print_log(1, f"Could not parse a molecule. (CID: {pc_id})")
            continue
        
        try:
            molecule_data: list = {
                "name": name,
                "id": pc_id,
                "mw": round(Descriptors.MolWt(molecule), 5),
                "h_donors": Descriptors.NumHDonors(molecule),
                "h_acceptors": Descriptors.NumHAcceptors(molecule),
                "log_p": round(Descriptors.MolLogP(molecule), 8)
            }

            violations, is_drug_likely = drug_candidate_check(molecule_data["mw"], molecule_data["h_donors"], molecule_data["h_acceptors"], molecule_data["log_p"])
            molecule_data.update({"is_drug_likely": is_drug_likely, "violations": violations})
            
            if record_full_data is True:
                molecule_data.update({"atom_number": molecule.GetNumAtoms(),
                                      "molecular_formula": rdMolDescriptors.CalcMolFormula(Chem.AddHs(molecule)),
                                      "bond_number": molecule.GetNumBonds(),
                                      "ring_number": molecule.GetRingInfo().NumRings(),
                                      "tpsa": Descriptors.TPSA(molecule),
                                      "rot_bond_number": Descriptors.NumRotatableBonds(molecule),
                                      "heavy_atoms": molecule.GetNumHeavyAtoms(),
                                      "aromatic_ring_number": Descriptors.NumAromaticRings(molecule),
                                      "sat_ring_number": Descriptors.NumSaturatedRings(molecule),
                                      "aliph_ring_number": Descriptors.NumAliphaticRings(molecule),
                                      "bridge_head_number": rdMolDescriptors.CalcNumBridgeheadAtoms(molecule),
                                      "sp3_percent": Descriptors.FractionCSP3(molecule),
                                      "surface_area": round(Descriptors.LabuteASA(molecule), 8),
                                      "partial_charge": round(Descriptors.MaxAbsPartialCharge(molecule), 8),
                                      "hk_alpha": round(Descriptors.HallKierAlpha(molecule), 8)
                                      })
            
            record_molecule_info(molecule_data)
            if save_csv_data is True:
                molecule_data["violations"] = "; ".join(molecule_data["violations"]) if molecule_data["violations"] else "None"
                append_csv_molecule_data(list(molecule_data.values()))
                
        except Exception as err:
            print_log(1, f"Compound descriptor error. {err}.")

    record_molecule_summary(id_list)
    if save_csv_data is True:
        save_csv_molecule_data()


def initialize_console_args() -> str:
    print("""\n\n\
 __       __            __                                __                           
/  \     /  |          /  |                              /  |                          
$$  \   /$$ |  ______  $$ |  ______    _______  __    __ $$ |        ______   __    __ 
$$$  \ /$$$ | /      \ $$ | /      \  /       |/  |  /  |$$ |       /      \ /  \  /  |
$$$$  /$$$$ |/$$$$$$  |$$ |/$$$$$$  |/$$$$$$$/ $$ |  $$ |$$ |      /$$$$$$  |$$  \/$$/ 
$$ $$ $$/$$ |$$ |  $$ |$$ |$$    $$ |$$ |      $$ |  $$ |$$ |      $$    $$ | $$  $$<  
$$ |$$$/ $$ |$$ \__$$ |$$ |$$$$$$$$/ $$ \_____ $$ \__$$ |$$ |_____ $$$$$$$$/  /$$$$  \ 
$$ | $/  $$ |$$    $$/ $$ |$$       |$$       |$$    $$/ $$       |$$       |/$$/ $$  |
$$/      $$/  $$$$$$/  $$/  $$$$$$$/  $$$$$$$/  $$$$$$/  $$$$$$$$/  $$$$$$$/ $$/   $$/ 
\n""")
    print(f"{'—'*90}\n  MolecuLex\n    Version: {application_version}\n    Made by Rim032\n    Last Updated: 2026/1/25\n{'—'*90}\n\n")
    
    console_parser = argparse.ArgumentParser()
    console_parser.add_argument("--fmin", help="Minimum value of CIDs to search through.")
    console_parser.add_argument("--fmax", help="Maximum value of CIDs to search through.")
    console_parser.add_argument("--file", help="The file path of a .txt containing CIDs formatted with commas or spaces.")
    console_parser.add_argument("--entry", help="A manually-entered string of CIDs.")

    console_parser.add_argument("--api_batch", help="Override the delay (in seconds) between API request batches. Be careful, this may ruin your work!")
    
    console_parser.add_argument("--noprint", action="store_true", help="Abstain from printing fetched data in console.")
    console_parser.add_argument("--nostat", action="store_true", help="Provide a short statistical summary regarding all of the fetched compounds.")
    
    console_parser.add_argument("--format", help="Display a short example of PubChem ID entry and file formatting for this program.")
    console_parser.add_argument("--full", action="store_true", help="Record all substansial properties associated with a compound.")
    console_parser.add_argument("--save_csv", help="Save the fetched data in a CSV file.")
    
    
    console_args = console_parser.parse_args()
    return console_args


def organize_args_input() -> str:
    user_input: str = ""
    
    if args.file is not None:
        user_input = read_id_file(args.file)
        user_input = format_id_list(user_input)
    if args.fmax is not None and args.fmin is not None:
        user_input = create_local_id_list(int(args.fmin), int(args.fmax))
    if args.entry is not None:
        user_input = format_id_list(args.entry)

    if args.api_batch is not None:
        global max_request_size
        max_request_size = int(args.api_batch)
        print_log(2, "Maximum API batch size has been changed. This may ruin your work!")

    if args.format:
        print("Format your file with each ID separated by a comma or space. (Ex: 256,438,1024)(Ex 2: 133 194 23 5)")
    if args.save_csv is not None:
        global save_csv_data
        save_csv_data = True

        global csv_file_name
        csv_file_name = args.save_csv
    if args.full:
        global csv_data
        csv_data = [["Compound Name", "PubChem ID", "Molecular Weight (g/mol)", "H-Donors", "H-Acceptors", "LogP", "Drug Possibility",
                     "Lipinski Violations", "Number of Atoms", "Molecular Formula", "Number of Bonds", "Number of Rings", "Topological Polar Surface Area",
                     "Rotatable Bonds", "Number of Heavy Atoms", "Number of Aromatic Rings", "Number of Saturated Rings", "Number of Aliphatic Rings",
                     "Number of Bridge Head Atoms", "Percent of sp^3 Hybridization", "Surface Area", "Max Partial Charge", "Hall Kier Alpha"]]

        global record_full_data
        record_full_data = True

    if args.nostat:
        global disable_stat_report
        disable_stat_report = True
    if args.noprint:
        global console_noprint
        console_noprint = True

        
    return user_input


def shutdown_cli() -> None:
    print_log(0, "No CID data was entered by the user... Type --help for more information.")
    print_log(0, "Exiting program...")

    time.sleep(4)
    sys.exit()



if __name__ == "__main__":
    application_start_time = time.perf_counter()
    
    args = initialize_console_args()
    user_input: str = organize_args_input()
    if user_input is None or user_input == "":
        shutdown_cli()

    smiles_list, compound_names_list = get_id_smiles(user_input)
    parse_compound_properties(smiles_list, compound_names_list, user_input)

    application_end_time = time.perf_counter()
    print_log(0, f"Search and analysis completed (Execution Time: {(application_end_time-application_start_time):2f} seconds)...")
    temp_char = input("\nPress enter to exit...")    

