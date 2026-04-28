  
**DOMAS: Domain Oriented Mapping of Alternative Splicing**

**Domain Oriented Mapping of Alternative Splicing (DOMAS), is a novel computational framework designed to bridge the gap between differential splicing events and their protein level effect. DOMAS accepts as input a list of differential splicing events, and performs a coordinate-aware mapping of these events onto the protein domain architectures. DOMAS then annotates each event with the affected domain(s), and classifies the effect as causing domain loss, gain, alteration, truncation or elongation, compared to the protein encoded by the canonical transcript.**

**Key Features**

* **Input formats:** Output of multiple differential splicing analysis, including LeafCutter (Li et al., 2016), rMATS2 (Wang et al., 2024), and MAJIQ (Vaquero-Garcia et al., 2023). To accommodate additional input formats, please write to us.  
* **Output format:** A csv list of domains affected by each alternative splicing event, and the class of the effect – loss, gain, alteration, elongation or truncation. Each event is linked to its visualization in DoCHAP format \[ref\].   
* **Protein domain databases:** DOMAS builds upon the DoChaP database (Gal-Oz et al., 2021), which includes domain annotation from list of databases and their reference.  
* Stay tuned \- DOMAS will also be available as a WebServer soon. 

**Installation & Requirements**

* Python 3.x  
* Dependencies: pandas, numpy, argparse, sqlite3, openpyxl (for Excel input) and matplotlib  
* **Database:** Requires access to a local instance of the **DoChaP DB**. See (Gal-Oz et al., 2021\) for installation instructions.

**Usage**

Run the utility from the command line:

Bash

python3 [domas.py](http://domas.py)  \-lc\_sig \<path-to-leafcutter\_ds\_significane.txt\> \-lc\_effect \<path-to-leafcutter\_ds\_effect\_sizes.txt file\> \-dochap \<dochap-path\> \-output\_csv \<output\_name.csv\>

**Input Parameters**

| Parameter | Description |
| :---- | :---- |
| \-lc\_sig | Path to leafcutter\_ds\_significane.txt output file |
| \-lc\_effect | Path to leafcutter\_ds\_effect\_sizes.txt output file |
| \-dochap | Path to the DoChaP DB directory/file. |
| \-output\_csv | The destination path for the generated results. |

**Output format**

The resulting CSV file provides:

* **comparison\_results:** Specific ID of the affected domain and the predicted change (e.g., gain/loss).  
* **domain\_descriptions:** Functional description of the affected domain retrieved from DoChaP.  
* Identification columns taken from the input**:** Includes cluster, gene\_name, gene\_id, start and end columns.

​References

​Gal-Oz, S. T., Haiat, N., Eliyahu, D., Shani, G., & Shay, T. (2021). DoChaP: the domain change presenter. *Nucleic Acids Research, 49*(W1), W162–W168. 

​Li, Y. I., Knowles, D. A., & Pritchard, J. K. (2016). LeafCutter: annotation-free quantification of RNA splicing. *Biorxiv,* , 044107\. 

​Vaquero-Garcia, J., Aicher, J. K., Jewell, S., Gazzara, M. R., Radens, C. M., Jha, A., Norton, S. S., Lahens, N. F., Grant, G. R., & Barash, Y. (2023). RNA splicing analysis using heterogeneous and large RNA-seq datasets. *Nature Communications, 14*(1), 1230\. 

​Wang, Y., Xie, Z., Kutschera, E., Adams, J. I., Kadash-Edmondson, K. E., & Xing, Y. (2024). rMATS-turbo: an efficient and flexible computational tool for alternative splicing analysis of large-scale RNA-seq data. *Nature Protocols, 19*(4), 1083–1104. 

​