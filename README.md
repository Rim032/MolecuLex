# MolecuLex (v0.9)
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
| `--drug` | Enables Lipinski’s Rule of Five screening and displays violations. |
| `--save_csv [name]` | Exports results to a CSV file. Provide a name or leave blank for default. |
| `--noprint` | Disables visual console output (useful for high-volume processing). |

### 2. API Flags
| Argument | Description |
| :--- | :--- |
| `--ovr_delay [int]` | Ovveride the delay between API request in seconds. Use with caution! |
| `--ovr_maxreq [int]` | Alter the maximum number of API requests before a cooldown engages. Use with caution! |
| `--ovr_cooldown [int]` | Change the temporary cooldown time after a limit of API requests is reached. Use with caution! |

---
**Example Usage & Arguments:**
Scan a range of IDs, evaluate drug-likeness, and export to a specific CSV file.
```bash
python moleculex.py --fmin 1 --fmax 512 --drug --save_csv my_data_and_such
```

---

## Misc.
* **Note:**
* **Testing:**
* **Compatibility:**
* **Notice:**
