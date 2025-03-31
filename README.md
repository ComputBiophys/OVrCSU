# OVrCSU

OVrCSU is a Python tool for calculating OVrCSU contacts in protein structures. This local implementation provides an alternative when the web server is unavailable.

## Installation

### Prerequisites
- Python 3.8 or later

### Dependencies
- MDAnalysis 
- numpy

## Usage
Basic command:
```bash
python OVrCSU.py -f input.pdb -o output.map
```

## Validation & Accuracy
Our implementation matches the web server results with <1% deviation. Comparative results with web server are provided below:

| Protein     | StateA               | StateB               | Contacts         | Multi_Contacts    |
|:-----------:|:--------------------:|:--------------------:|:----------------:|:-----------------:|
| GlnBP       | 505/505              | 552/552              | 440/440          | 177/177           |
| Arc         | 190/190              | 175/176 (+1)         | 96/96            | 173/174 (+1)      |
| Arc (Other) | 190/184 (-8/+2)      | 175/166 (-11/+2)     | 96/91 (-5)       | 173/168 (-8/+13)  |
| Hinge       | 346/347 (+1)         | 372/374 (+2)         | 254/254          | 210/213 (+3)      |
| SemiSWEET   | 336/336              | 344/344              | 252/252          | 176/176           |
| TRAAK       | 1026/1030 (+4)       | 1035/1036 (+1)       | 918/921 (+3)     | 225/224 (-2/+1)   |

**Note:** 
- StateA and StateB refer to contacts in the single-basin Go-Martini model with proteins as State A and B (cutoff: 0.3–1.2 nm, 3-residue gap exclusion).  
- Contacts and Multi_Contacts refer to shared contacts and state-unique contacts in the multiple-basin Go-Martini model, respectively.  
- Values are shown as Web Server/This Script. The number in the brackets indicates the difference between the web server and this script.
- "Other" refers to [GoMartini's ContactMapGenerator](https://github.com/Martini-Force-Field-Initiative/GoMartini/tree/main/ContactMapGenerator).


## References
Please cite these works if using this tool:  
[1] Wołek, K.; Gómez-Sicilia, À.; Cieplak, M. Determination of contact maps in proteins: A combination of structural and chemical approaches. The Journal of Chemical Physics 2015, 143.  
[2] Yang, S.; Song, C. Multiple-basin Go-Martini for investigating conformational transitions and environmental interactions of proteins. bioRxiv 2024.  
