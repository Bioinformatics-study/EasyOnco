### 🎨 EasyOnco

## 💬 Introduction
EasyOnco is a tool that create an MAF file with the result file of [PiSeq](https://www.annlabmed.org/journal/view.html?doi=10.3343/alm.2023.43.4.328) pipeline and draw oncoplot using maftools.

## ⚙️ Requirements
Python version >= 3.0

Python libraries
- [Pandas](https://pypi.org/project/pandas/)
- [Argparse](https://pypi.org/project/argparse/)
- [openpyxl](https://openpyxl.readthedocs.io/en/stable/tutorial.html)
```
# install pandas
pip install pandas

# install argparse
pip install argparse

# install openpyxl
pip install openpyxl
```

R library
- [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)(version > 2.20.0)
- [readxl](https://cran.r-project.org/web/packages/readxl/readme/README.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/readme/README.html)
- [scales](https://cran.r-project.org/web/packages/scales/readme/README.html)
- [circlize](https://github.com/jokergoo/circlize)

```
# install maftools
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("maftools")

# install readxl
install.packages("readxl")

# install dplyr
install.packages("dplyr")

# install scales
install.packages("scales")

# install circlize
install.packages("circlize")
```

## ⬇️ Installation
```
git clone https://github.com/Bioinformatics-study/EasyOnco.git
```

## 💻 Usage
> EasyOnco
```
python EasyOnco.py [-i] [-f] [-s] [-o]
```
```
python EasyOnco.py -i *xlsx -f v vv x -s P -o example.maf
```

*Options*
- `-f` : Filtering field, enter the character marked as true.
- `-s` : Fill in the Excel sheet after conducting mutation analysis. [P: path.vus, A: All variants]
- `-i` : Specify the xlsx file for each sample to create the Oncoplot (Use regular expressions if possible).
- `-o` : Specifies the name of the maf file to be generated as the result. If not specified, the result is generated under the name EasyOnco.maf by default.

## 📊 Output
### Oncoprint & mafsummary plot
![그림1](https://github.com/user-attachments/assets/140e8f04-2526-4a14-a216-5b4798254f49)
### lollipop & vaf plot
![2](https://github.com/user-attachments/assets/a6a3e34b-b1fa-4340-9986-a4de4ce9e41b)
### ti/tv & druggable plot
![3](https://github.com/user-attachments/assets/b6b18976-a730-4dd0-8c1b-c2d90c3bdb26)

## 🤖 Running Manually
> Steps to create the default output
1. Set the options in the `oncoplot_options_v1.xlsx` file according to the analysis requirements
2. Navigate to the directory containing the xlsx file with mutation information.
3. If annotation information is available, save it as a file named `Clinical_annotation.txt`.
4. Enter the Input file, filter string, and xlsx sheet information as shown in the Python code example above.

## 🌈 For Customization
> Steps to customize colors and annotation bars, etc
1. Set the options in the `oncoplot_options_v1.xlsx` file.
2. Sections written in blue text can be customized with values of your choice.
3. Black text indicates that only values set in the drop down can be selected.
4. Clicking each cell will display a description of the values that can be entered.
5. You can refer to the option sheet for the color code.

## 👥 Authors
[YUHS Labmed](https://sites.google.com/view/diagnosticlaboratory/home)
