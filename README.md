# CDAtransformer

**CDAtransformer** is an R Shiny application for parsing **C-CDA (CDA XML)** and **FHIR (JSON)** clinical documents into analysis-ready tabular outputs.

The application supports both single-document uploads and table-based datasets where each row contains a complete document payload associated with a patient identifier.

![CCDASample](https://github.com/BoyceLab/CDAtransformer/blob/main/CCDASample.jpg)

---

## Overview

Clinical data are frequently distributed as nested XML or JSON documents (C-CDA or FHIR), which are difficult to analyze directly in R. `CDAtransformer` extracts structured elements from these documents and returns tidy tables suitable for filtering, joins, and downstream analysis.

This tool is intended for research, informatics, and data science workflows using exported electronic health record (EHR) data.

---

## Features

- Parse **C-CDA (CDA XML)** documents
- Parse **FHIR JSON**, including Bundle and `searchset` responses
- Accept either:
  - Single document uploads, or
  - CSV datasets containing one full document per row
- Preserve patient identifiers in all outputs
- Multiple export formats:
  - Long (tidy, recommended)
  - Wide (single-row summary)
  - FHIR split export (one CSV per resource type)

---

## Supported Input Formats

### C-CDA (CDA XML)

#### Single document upload

Upload a `.txt` file containing **one complete C-CDA XML document**.

#### Dataset upload (CSV)

The CSV file must contain **exactly two columns**:

| Column name   | Description |
|-------------|-------------|
| `patient_id` | Identifier used to link results back to the patient |
| `doc`        | Full C-CDA XML document as a single string |

**Important:**  
Each row in `doc` must contain **one complete XML document**. If XML content is split across rows due to quoting or newline issues, parsing will fail.

---

### FHIR

#### Single document upload

Upload a `.txt` file where **each non-empty line** is a valid JSON object:
- a FHIR resource, or
- a FHIR Bundle (e.g., `searchset`)

#### Dataset upload (CSV)

The CSV file must contain **exactly two columns**:

| Column name    | Description |
|--------------|-------------|
| `patient_id`  | Identifier used to link results back to the patient |
| `record_body` | JSON string containing a FHIR resource or Bundle |

FHIR Bundles with `entry[].resource` elements are supported.

---

## Output Formats

### Long format (recommended)

A tidy, row-based output suitable for analysis.

### Wide format

Collapses repeated paths into a single row per document. Intended for inspection only; not recommended for analysis.

### FHIR split export

Exports one CSV file per FHIR resource type (e.g., `Patient.csv`, `Condition.csv`, `Observation.csv`), packaged as a ZIP archive.

Each file includes `patient_id` for linkage.

---

## File Size Limits

Upload size is controlled in `app.R`:

```r
options(shiny.maxRequestSize = 250 * 1024^2)  # 250 MB
```

When deployed behind a proxy (e.g., nginx, Shiny Server, Posit Connect), additional upload limits may apply outside of R.

---

## Installation

### Prerequisites

- R (>= 4.0)
- RStudio (recommended)

### Install dependencies

<details>
<summary>Show R commands</summary>

```r
install.packages(c(
  "shiny", "xml2", "jsonlite", "dplyr", "tidyr", "DT",
  "tibble", "zip", "readr", "tools"
))
```
</details>

### Install from GitHub

<details>
<summary>Show R commands</summary>

```r
install.packages("devtools")
devtools::install_github("BoyceLab/CDAtransformer", force = TRUE)
```
</details>

---

## Usage

Run the Shiny application:

<details>
<summary>Show R commands</summary>

```r
library(shiny)
runApp(".")
```
</details>

---

## Typical dataset workflow

1. Export clinical documents to CSV
2. Confirm required columns:
   - **C-CDA:** `patient_id`, `doc`
   - **FHIR:** `patient_id`, `record_body`
3. Upload the CSV into the application
4. Select document type (**C-CDA** or **FHIR**)
5. Choose an export format
6. Download the results

---

## Contributing

Contributions are welcome. Please open an issue or submit a pull request via GitHub.

---

## License

MIT License. See the `LICENSE` file for details.

