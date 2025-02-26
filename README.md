# shinyGenePlot
A simple Shiny application for visualizing gene expression data across different experimental conditions.

## Online version
https://molinerislab.shinyapps.io/genePlots/

## Features

- Upload and visualize gene expression data with an intuitive interface
- Select specific genes of interest for visualization
- Group samples by any metadata variable (e.g., condition, treatment, time point)
- Reorder factor levels via drag-and-drop interface
- Color boxplots by any metadata variable
- View data in full-screen plot mode
- Inspect data through tabular previews

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/gene-expression-visualization.git
cd gene-expression-visualization
```

2. Make sure you have R installed, then install the required packages:
```r
install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "readr", "DT", "sortable"))
```

3. Run the app:
```r
# From R console
shiny::runApp()

# Or from terminal
Rscript -e "shiny::runApp()"
```

## Input Data Format

The app requires two input files:

### 1. Gene Expression Matrix
- CSV, TSV, or gzipped (.gz) file
- Genes in rows, samples in columns
- First column contains gene identifiers
- Example:

```
gene,sample1,sample2,sample3
GENE1,1.2,2.3,3.4
GENE2,2.2,3.3,4.4
GENE3,3.2,4.3,5.4
```

### 2. Sample Metadata
- CSV, TSV, or text file
- Rows represent samples
- First column contains sample identifiers that match expression matrix column names
- Remaining columns contain metadata (condition, treatment, etc.)
- Example:

```
sample_id,condition,cell_type
sample1,treated,neuron
sample2,control,neuron
sample3,treated,glia
```

## Usage

1. Upload your gene expression matrix and sample metadata files
2. Select genes of interest from the dropdown
3. Choose which metadata column to use for grouping samples
4. Drag and reorder factor levels to customize the plot
5. Select which metadata column to use for coloring
6. Click "Update Plot" to generate the visualization
7. Switch between tabs to view the plot or inspect data

## Example

When correctly set up, the app will generate boxplots for each selected gene, grouped by the chosen metadata variable, and colored by another variable of your choice.

## Notes and Troubleshooting

- Sample IDs must match exactly between your expression matrix and metadata
- If some samples appear as "NA" in the plot, check for mismatches between files
- The app automatically trims whitespace from sample IDs to prevent join issues
- You will receive a notification if samples in your expression data are missing from metadata
- For large datasets, the initial loading may take a few moments

## License

AGPL3

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
