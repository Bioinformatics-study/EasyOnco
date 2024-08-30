## 🎨 EasyOnco

## 💬 Introduction
EasyOnco is a tool that create an MAF file with the result file of [PiSeq](https://www.annlabmed.org/journal/view.html?doi=10.3343/alm.2023.43.4.328) pipeline and draw oncoplot using maftools.



## 👥 Authors
[YUHS Labmed](https://sites.google.com/view/diagnosticlaboratory/home)

## ⚙️ Requirements
Python version >= 3.0

Python libraries
- [Pandas](https://pypi.org/project/pandas/)
- [Argparse](https://pypi.org/project/argparse/)

R library
- [maftools](https://bioconductor.org/packages/release/bioc/html/maftools.html)
- readxl
- dplyr
- scales
- circlize

```
# install pandas
pip install pandas

# install argparse
pip install argparse

# install maftools
R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("maftools")
```

## ⬇️ Installation
```
git clone https://github.com/Bioinformatics-study/EasyOnco.git
```

## 💻 Usage
> OnGoMAF
```
python OnGoMAF.py [-f] [-s] [-i]
```
```
python EasyOnco.py -i *xlsx -f v vv x -s P
```

*Options*
- `-f` : Filtering field, enter the character marked as true.
- `-s` : Fill in the Excel sheet after conducting mutation analysis. [P: path.vus, A: All variants]
- `-i` : Specify the xlsx file for each sample to create the Oncoplot (Use regular expressions if possible).

## 📊 Output
1. Oncoprint
![240828_oncoplot](https://github.com/user-attachments/assets/ad6a96ad-8d30-4d7f-a567-d9a6c8692038)

2. Variant summary plot
![240828_mafsummary](https://github.com/user-attachments/assets/0d1baed5-14cb-48b3-8d1d-3e4922576325)

3. ti/tv plot
![240828_onco titv](https://github.com/user-attachments/assets/2a2ef071-ad9d-42b3-9d18-bf602cd8e55e)

4. vaf boxplot
![240828_vafplot](https://github.com/user-attachments/assets/f3f44f7a-9993-4dd1-84a4-2e5c2794260b)

5. lollipop plot
![240828_lollipop](https://github.com/user-attachments/assets/4bad8f88-d888-41c3-97a2-2ce17c911bd6)

6. druginteraction bartplot
![240828_druginteraction](https://github.com/user-attachments/assets/7871e047-9ade-4ffe-ae9f-47de2735529d)

## 🤖 Running Manually
> Step to create the default output
1. Set the options in the `oncoplot_options_v1.xlsx` file according to the analysis requirements
2. Navigate to the directory containing the xlsx file with mutation information.
3. If annotation information is available, save it as a file named `Clinical_annotation.txt`.
4. Enter the Input file, filter string, and xlsx sheet information as shown in the Python code example above.

## 🌈 For Customization
> Step to 
