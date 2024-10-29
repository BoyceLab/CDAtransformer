# Boyce Lab C-CDA and FHIR Data Transformer

This Shiny application allows users to upload and parse C-CDA and FHIR text files. It extracts relevant elements and presents them in a tabular format for easy viewing and download. 

![Screenshot of the Boyce Lab C-CDA and FHIR Data Converter](https://github.com/BoyceLab/CDAtransformer/blob/main/Screenshot%202024-06-23%20203614.jpg)

## Features

- Upload C-CDA and FHIR text files (up to 30MB)
- Parse and display the contents of C-CDA and FHIR files
- Download the parsed data as a CSV file

## Installation

### Prerequisites

Make sure you have R and RStudio installed on your system. You will also need the `devtools` package to install the app from GitHub.

1. **Install R:** [Download R](https://cran.r-project.org/)
2. **Install RStudio:** [Download RStudio](https://www.rstudio.com/products/rstudio/download/)

### Install Dependencies

Open RStudio and run the following commands to install the necessary packages:

```r
install.packages(c("shiny", "xml2", "jsonlite", "dplyr", "tidyr", "DT", "devtools"))
```
### Clone the Repository from GitHub
Open a terminal or command prompt and navigate to the directory where you want to clone the repository. Then, run the following command:
```r
git clone https://github.com/BoyceLab/CDAtransformer.git
```
### Install the App from GitHub

Use the `devtools` package to install the app from GitHub:

```r
devtools::install_github("BoyceLab/CDAtransformer", force = TRUE)
```

## Usage

1. **Set the Working Directory:**

   Open RStudio and set the working directory to where your `app.R` file is located:

   ```r
   setwd("path/to/CDAtransformer")
   ```

   Replace `"path/to/CDAtransformer"` with the actual path to the cloned repository.

2. **Run the App:**

   Load the required packages and run the Shiny app:

   ```r
   library(shiny)
   runApp(".")
   ```

3. **Upload a File:**

   - Click the "Choose CDA or FHIR Text File" button to upload your text file.
   - Select the file type (CDA or FHIR).
   - Click the "Parse Document" button to parse and display the contents.

4. **Download Parsed Data:**

   - Click the "Download CSV" button to download the parsed data as a CSV file.

## Code Overview

### UI

The UI is defined in the `ui` function and includes a title panel, file input, radio buttons for file type selection, a parse button, and a download button. The parsed data is displayed in a datatable.

### Server

The server logic is defined in the `server` function and includes:

- **File Parsing:** Functions to parse CDA and FHIR files.
- **Data Display:** Render the parsed data in a datatable.
- **Data Download:** Allow users to download the parsed data as a CSV file.

### Helper Functions

- `parse_cda(file)`: Parses a CDA file and extracts relevant elements.
- `flatten_list(x, name)`: Recursively flattens a list.
- `parse_fhir(file)`: Parses multiple FHIR resources and flattens them.
- `transpose_data(data)`: Transposes a data frame.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

