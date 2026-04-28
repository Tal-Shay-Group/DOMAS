**DOMAS: Domain Oriented Mapping of Alternative Splicing**

**DOMAS** is a computational framework designed to bridge the gap
between alternative splicing products and proteins. It maps alternative
splicing (AS) events onto protein domain architectures using a
coordinate-aware approach.

**Key Features**

**Functional Annotation:** Identifies domain gain, loss, alteration in sequence or length
caused by alternative splicing.

**Database Integration:** Utilizes the DoChaP database for high-fidelity
transcript and domain mapping.

**Installation & Requirements**

Python 3.x

Dependencies: pandas, numpy, argparse, sqlite3, openpyxl (for Excel
input) and matplotlib

**Database:** Requires access to a local instance of the **DoChaP DB**.

**Usage**

Run the utility from the command line:

Bash

python3 domas.py -input \<input_file.xlsx\> -dochap \<dochap-path\>
-output_csv \<output_name.csv\>

**Input Parameters**

  -------------------- --------------------------------------------------
  **Parameter**        **Description**

  -input               Path to an Excel file (.xlsx) containing junction
                       data.

  -dochap              Path to the DoChaP DB directory/file.

  -output_csv          The destination path for the generated results.
  -------------------- --------------------------------------------------

**Input File Format (Excel)**

The input file must contain the following columns:

**h_junction:** Junction coordinates in format Chr:Start:End (e.g.,
Chr6:33416776:33417096).

**ensembl_h:** Ensembl Gene ID (e.g., ENSG00000112514).

**cluster:** A unique identifier used to group junctions for comparison.

**rank_h:** The exons connected by the junction taken from the canonical
transcript. Format: \<start-exon\>\_\<end-exon\> (e.g. E2_E4.).Included
for readability; not used in the core analysis.

**symbol_h:** Gene common name. (e.g. CUTA). Included for readability;
not used in the core analysis.

**Output Columns**

The resulting CSV file provides:

**comparison_results:** Specific ID of the affected domain and the
predicted change (e.g., gain/loss).

**domain_descriptions:** Functional description of the affected domain
retrieved from DoChaP.

Identification columns taken from the input**:** Includes cluster,
gene_name, gene_id, start and end columns.
