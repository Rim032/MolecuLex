# MolecuLex (v0.99)
```
 __       __            __                                __                           
/  \     /  |          /  |                              /  |                          
$$  \   /$$ |  ______  $$ |  ______    _______  __    __ $$ |        ______   __    __ 
$$$  \ /$$$ | /      \ $$ | /      \  /       |/  |  /  |$$ |       /      \ /  \  /  |
$$$$  /$$$$ |/$$$$$$  |$$ |/$$$$$$  |/$$$$$$$/ $$ |  $$ |$$ |      /$$$$$$  |$$  \/$$/ 
$$ $$ $$/$$ |$$ |  $$ |$$ |$$    $$ |$$ |      $$ |  $$ |$$ |      $$    $$ | $$  $$<  
$$ |$$$/ $$ |$$ \__$$ |$$ |$$$$$$$$/ $$ \_____ $$ \__$$ |$$ |_____ $$$$$$$$/  /$$$$  \ 
$$ | $/  $$ |$$    $$/ $$ |$$       |$$       |$$    $$/ $$       |$$       |/$$/ $$  |
$$/      $$/  $$$$$$/  $$/  $$$$$$$/  $$$$$$$/  $$$$$$/  $$$$$$$$/  $$$$$$$/ $$/   $$/
```

**MolecuLex** is a versatile, light-weight chemoinformatics CLI tool designed for automated reagent, drug, and general molecule probing. Through the usage of **RDKit** and **PubChemPy**, researchers are able to rapidly evaluate large batches of compounds en masse.

![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Version](https://img.shields.io/badge/version-0.99-orange)

## Description

This program utilizes a high-volume pipeline optimized for chunked [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) requests at increments of 256 compounds. It hosts a robust chemoinformatics suite through providing metrics on molecular composition, topological polar surface area (TPSA), electronic partial charges, and so much more, alongside Lipinski Rule of Five screening. Data is able to be rapidly organized and exported through a .CSV file pipeline for molecule collection, database screening, and more. Supporting sequential ranges, manual entry, and file parsing, it serves as a flexible and easy-to-use tool for chemists and others alike.

## Commands & Settings

### 1. Data Input Flags (Choose one)
| Argument | Description | Usage Example |
| :--- | :--- | :--- |
| `--fmin [int]` & `--fmax [int]` | Scans a sequential range of PubChem CIDs. | `--fmin 100 --fmax 200` |
| `--file [name/path]` | Reads CIDs from a `.txt` file (comma or space separated). | `--file list.txt` |
| `--entry [CIDs]` | Allows manual entry of CIDs directly in the console. | `--entry 2244, 1983` |

### 2. Analysis & Export Flags
| Argument | Description |
| :--- | :--- |
| `--save_csv [name]` | Exports results to a CSV file. Provide a name or leave blank for default. |
| `--full` | Gathers additional data and properties for each compound. |
| `--noprint` | Disables console printing of each compound's properties (useful for high-volume processing). |
| `--nostat` | Disables console printing of a cumulative data summary. |
| `--api_batch [int]` | Override the amount of requests made in one chunk at a time. Use with caution! |

---
**Example Usage & Arguments:**
Scan a range of IDs, search for additional structural properties, and export to a specific CSV file.
```bash
python moleculex.py --fmin 1 --fmax 4096 --full --save_csv my_data_and_such
```

---

## Misc.
> [!NOTE]
> Files containing CIDs must be separated by commas (exclude any and all spaces. **type --format for more information in the program**).
* **Testing:** During testing, this program has been able to successfully analyze and export chemical data for over 120,000+ compounds and counting. Over eight trials, consisting of a thousand entries each, it has demonstrated to be capable of requesting and processing 130 (± 2.2) compounds per second on average.
* **Compatibility:** Windows 10/11 and [Python 3.10+](https://www.python.org/downloads/) are required.
