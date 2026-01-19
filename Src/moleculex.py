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


os.system('')
#PUG REST throttling to avoid request denial and server strain.
server_request_delay = 0.1
max_pc_requests = 200
max_request_cooldown = 30
desired_extension = ".txt"

log_types = [
    "LOG",      #0
    "ERROR",    #1
    "WARNING"   #2
]

display_drug_attributes = False
application_version = 0.9
console_noprint = False
save_csv_data = False


csv_file_name = "untitled_molecuLex_data"
csv_data = [
    ["Compound Name", "Molecular Weight (g/mol)", "H-Donors", "H-Acceptors", "LogP", "Drug Possibility", "Lipinski Violations", "PubChem ID"]
]




def print_log(log_type: int, message: str) -> str:
    if message is None or log_type is None:
        return
    print(f"[{log_types[log_type]}] - {message}")


def append_csv_molecule_data(molecule_data: list[any]) -> None:
    if molecule_data is None:
        print_log(1, f"Cannot append molecule data to CSV data list. Passed data is invalid.")
        return

    global csv_data
    csv_data.append(molecule_data)


def save_csv_molecule_data() -> None:
    if save_csv_data is False:
        print_log(1, f"Cannot save molecule CSV file.")
        return

    try:
        with open(f"{csv_file_name}.csv", mode='w', newline='', encoding='utf-8') as csv_file:
            file_writer = csv.writer(csv_file, delimiter=',')
            file_writer.writerows(csv_data)
        print_log(0, f"Successfully saved data to {csv_file_name}.csv")
    except Exception as err:
        print_log(1, f"Failed to save CSV: {err}")


def record_molecule_info(name: str, mw: float, hdonors: int, hacceptors: int, logp: float, drug_chance: bool, violations: list[str], pc_id: int) -> str:
    if console_noprint is True:
        return
    
    WIDTH = 62
    LABEL_W = 20
    DATA_W = 30
    
    try:
        symbol = "•" if drug_chance else "♦"
        tag = "[ LIKELY DRUG CANDIDATE ]" if drug_chance else "[ UNLIKELY DRUG CANDIDATE ]"
        print(f"\n╔{'═' * (WIDTH-2)}╗")
        
        prefix = f"{symbol} CID ({pc_id}): "
        available_space = (WIDTH - 4) - len(prefix)

        if len(name) > available_space:
            display_name = name[:available_space - 2] + ".."
        else:
            display_name = name
        full_content = f"{prefix}{display_name}"
        
        print(f"║ {full_content:<{WIDTH - 4}} ║")
        print(f"╠{'═' * (WIDTH-2)}╣")

        print(f"║ { 'Mass:':<{LABEL_W}}{mw:>{DATA_W}.4f} g/mol {' ':<2}║")
        print(f"║ { 'LogP:':<{LABEL_W}}{logp:>{DATA_W}.4f} {' ':<8}║")
        print(f"║ { 'H-Donors:':<{LABEL_W}}{hdonors:>{DATA_W}} {' ':<8}║")
        print(f"║ { 'H-Acceptors:':<{LABEL_W}}{hacceptors:>{DATA_W}} {' ':<8}║")

        if display_drug_attributes is True:
            print(f"╟{'─' * (WIDTH-2)}╢")
            print(f"║ {tag:^{WIDTH-4}} ║")
            
            if violations:
                for v in violations:
                    clean_v = (v[:WIDTH-10] + '..') if len(v) > (WIDTH-10) else v
                    print(f"║ ! {clean_v:<{WIDTH-5}}║")
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


def get_id_smiles(id_list: list[str]) -> list[str]:
    if id_list is None:
        print_log(1, "An invalid list of formatted IDs was passed.")
        return

    smiles_list: list[str] = []
    name_list: list[str] = []

    bar_length = 30
    request_counter = 0
    try:
        print(f"\n{'—'*20} DATA ACQUISITION PHASE {'—'*20}")
        print_log(0, f"Connecting to PubChem API for {len(id_list)} compounds...")
        for i, pub_id in enumerate(id_list, 1):
            if request_counter > max_pc_requests:
                print("")
                print_log(2, f"Currently cannot faciliate any additional request to PubChem. ({max_pc_requests} request MAX)")
                print_log(0, f"Requests have temporarily been placed on cooldown for {max_request_cooldown} seconds.")
                
                time.sleep(max_request_cooldown)
                request_counter = 0
                    
            pub_compound = pcp.Compound.from_cid(pub_id)
            time.sleep(server_request_delay)

            progress = i / len(id_list)
            filled_len = int(bar_length * progress)
            bar = "▓" * filled_len + "░" * (bar_length - filled_len)
            percent = int(progress * 100)
            sys.stdout.write(f"\rPROGRESS: [{bar}] {percent}% | Analyzing CID: {pub_id:<10}")
            sys.stdout.flush()
            
            if pub_compound is not None:
                smiles_list.append(pub_compound.smiles)
                name_list.append(pub_compound.iupac_name)
            request_counter = request_counter + 1
    except Exception as err:
        print("")
        print_log(1, f"API Issue. {err}.")
    print("")

    return smiles_list, name_list


def check_drug_properties(smiles_list: list[str], names_list: list[str], id_list: list[str]) -> None:
    if smiles_list is None:
        print_log(1, "An invalid SMILES list was passed.")
        return
    
    for smiles, name, pc_id in zip(smiles_list, names_list, id_list):
        molecule = Chem.MolFromSmiles(smiles)

        try:
            mw: float = round(Descriptors.MolWt(molecule), 5)
            h_donors: int = Descriptors.NumHDonors(molecule)
            h_acceptors: int = Descriptors.NumHAcceptors(molecule)
            logp: float = round(Descriptors.MolLogP(molecule), 8)
            is_drug_likely = True


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
                is_drug_likely = False

            record_molecule_info(name, mw, h_donors, h_acceptors, logp, is_drug_likely, violations, pc_id)
            if save_csv_data is True:
                violations_formatted = "; ".join(violations) if violations else "None"
                append_csv_molecule_data([name, mw, h_donors, h_acceptors, logp, is_drug_likely, violations_formatted, pc_id])
        except Exception as err:
            print_log(1, f"Compound descriptor error. {err}.")

    if save_csv_data is True:
        save_csv_molecule_data()


def initialize_console_args() -> str:
    print("""\n\n
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
    print(f"{'—'*90}\n  MolecuLex\n    Version: {application_version}\n    Made by Rim032\n    Last Updated: 2026/1/18\n{'—'*90}\n\n")
    
    console_parser = argparse.ArgumentParser()
    console_parser.add_argument("--fmin", help="Minimum value of CIDs to search through.")
    console_parser.add_argument("--fmax", help="Maximum value of CIDs to search through.")
    console_parser.add_argument("--file", help="The file path of a .txt containing CIDs formatted with commas or spaces.")
    console_parser.add_argument("--entry", help="A manually-entered string of CIDs.")

    console_parser.add_argument("--ovr_delay", help="Override the delay time between PubChem requests. Be careful, this may ruin your work!")
    console_parser.add_argument("--ovr_maxreq", help="Override the number of maximum requests in an operation period. Be careful, this may ruin your work!")
    console_parser.add_argument("--ovr_cooldown", help="Override the time of a temporary cooldown period. Be careful, this may ruin your work!")
    
    console_parser.add_argument("--drug", action="store_true", help="Record and print drug viability and violations.")
    console_parser.add_argument("--noprint", action="store_true", help="Abstain from printing fetched data in console.")
    console_parser.add_argument("--format")
    
    console_parser.add_argument("--full", help="Record all substansial properties associated with a compound.")
    console_parser.add_argument("--csv_save", help="Save the fetched data in a CSV file.")
    
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

    if args.ovr_delay is not None:
        global server_request_delay
        server_request_delay = float(args.ovr_delay)
        print_log(2, "Maximum request delay has been changed. This may ruin your work!")
    if args.ovr_maxreq is not None:
        global max_pc_requests
        max_pc_requests = int(args.ovr_maxreq)
        print_log(2, "Maximum period requests has been changed. This may ruin your work!")
    if args.ovr_cooldown is not None:
        global max_request_cooldown
        max_request_cooldown = float(args.ovr_cooldown)
        print_log(2, "Maximum request cooldown time has been changed. This may ruin your work!")

    if args.drug:
        global display_drug_attributes
        display_drug_attributes = True
    if args.format:
        print("Format your file with each ID separated by a comma or space. (Ex: 256,438,1024)(Ex 2: 133 194 23 5)")
    if args.noprint:
        global console_noprint
        console_noprint = True
    if args.csv_save is not None:
        global save_csv_data
        save_csv_data = True

        global csv_file_name
        csv_file_name = args.csv_save
        

    return user_input



if __name__ == "__main__":
    application_start_time = time.perf_counter()
    
    args = initialize_console_args()
    user_input: str = organize_args_input()
    
    if user_input is None or user_input == "":
        print_log(0, "No CID data was entered by the user... Type help for more information.")
        print_log(0, "Exiting program...")

        time.sleep(3)
        sys.exit()
        
    compound_lists: list[str] = get_id_smiles(user_input)
    check_drug_properties(compound_lists[0], compound_lists[1], user_input)

    application_end_time = time.perf_counter()
    print_log(0, f"Search and analysis completed (Execution Time: {(application_end_time-application_start_time):2f} s)...")
    temp_char = input("\nPress enter to exit...")    
