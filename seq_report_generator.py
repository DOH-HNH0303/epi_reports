from fpdf import FPDF
from PyPDF2 import PdfFileReader, PdfReader, PdfWriter
import numpy as np
import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image


def add_dendrogram_as_pdf(pdf_name, tree_file, format="newick", output_filename=""):
  tree = Phylo.read(tree_file, format)
  tree.ladderize()
  count = 0
  for i in tree.get_terminals():
    count+=1
  
  mdist = max([tree.distance(tree.root, x) for x in tree.get_terminals()])
  plt.rc('font', size=7)
  plt.rc('lines', linewidth=0.5)
  fig = plt.figure()
  fig.set_size_inches(10, count*0.21)
  axes = fig.add_subplot(1, 1, 1)
  
  Phylo.draw(tree, 
             axes=axes,#ax, 
             do_show=False, 
             show_confidence=False,
             xticks=([],), yticks=([],),
             ylabel=('',), xlabel=('',),
             xlim=(-mdist*0.1,mdist+mdist*0.1),
             axis=('off',),)
  if output_filename:
    base=output_filename
  else:
    base=(Path(tree_file).stem)+".pdf"
  plt.savefig(base, dpi=150, bbox_inches='tight')
  
  #add_existing_pdf(pdf_name, base)

#def add_existing_pdf(pdf_name, existing_file):
 # pdf_name.Template(existing_file, x = 10, y = 10, w = 80, h = 30, type = 'P')
def add_image(pdf_name, file):
  pdf_name.ln()
  img = Image.open(file)
  width, height = img.size
  new_height = round((200 * height) / width)
  pdf_name.image("pangenome_matrix.png", h=new_height, w=200)
  

def add_page_header(pdf_name, text='Header of PDF Report', font='Times', fontsize=16, bold=True,underline=False,italic=False):
  """Creates new page in pdf and adds text as header"""
  pdf_name.add_page()
  font_style = _font_style(bold, underline, italic)
  pdf_name.set_font(font, font_style, fontsize)
  w = pdf_name.get_string_width(text) + 6
  pdf_name.set_x((210 - w) / 2)
  pdf_name.cell(w, 9, text, 0, 0, 'C')
  pdf_name.line(20, 18, 210 - 20, 18)

def add_section_header(pdf_name,text = "Section Header", font='Times', fontsize=12, bold=True,underline=False,italic=False):
  pdf_name.ln()
  font_style = _font_style(bold, underline, italic)
  pdf_name.set_font(font, font_style, fontsize)
  pdf_name.set_fill_color(200, 220, 255)
  pdf_name.cell(0, 6, text, 0, 1, 'L', 1)

def add_paragraph(pdf_name, text="", font="Times", fontsize=12, bold=False,underline=False,italic=False):
  pdf_name.ln(10)
  font_style = _font_style(bold, underline, italic)
  pdf_name.set_font(font, font_style, fontsize)
  pdf_name.multi_cell(0, 5, text)

def add_table(pdf_name, df, col_len_dict={0:24, 1:20}, max_line_len = 40):
  pdf_name.ln(10)
  pdf_name.set_font("Times", size=10)
  line_multiplier = pdf_name.font_size * 2.3
  line_height = pdf_name.font_size * 2.5
  data, width_ratios, max_height_list = _reformat(df, max_line_len=max_line_len)
  df_cp = df.copy()
  col_width = []
  for i in width_ratios:
    col_width.append(pdf_name.epw*i/sum(width_ratios))

  hold = col_width
  lh_list = [] #list with proper line_height for each row

  count = 0 
  for index, row in df_cp.iterrows():
      new_line_height = 7
      for datum in row:
          datum = str(datum)
          word_list = datum.split()
          number_of_words = len(word_list) #how many words
          if number_of_words>2 or len(datum)>10:#names and cities formed by 2 words like Los Angeles are ok)
              l_height = pdf_name.font_size * (number_of_words/2)+1 #new height change according to data 
              if l_height > new_line_height:
                new_line_height = l_height  #new height change according to data 
      count +=1
      lh_list.append(new_line_height)

  #create your fpdf table ..passing also max_line_height!
  if col_len_dict:      
    for key in list(col_len_dict.keys()):
      hold[key] = col_len_dict[key]
  hold[0]=24
  hold[1]=20

  pdf_name.set_fill_color(r=226, g=228, b=249)
  # Create table headers
  for idx in range(len(df.columns)):
    width = hold[idx]
    col = df.columns[idx]
    line_height = line_multiplier
    pdf_name.multi_cell(width, line_height, col, border=1,align='L',
                  new_x="RIGHT", new_y="TOP", max_line_height=pdf_name.font_size, fill=True)
  pdf_name.ln(line_height)

  # Add data to table
  for j,row in enumerate(data):
      if j != len(lh_list):
        for idx in range(len(row)):
            datum = str(row[idx])
            width = hold[idx]
            line_height = lh_list[j] #choose right height for current row
            pdf_name.multi_cell(width, line_height, datum, border=1,align='L',ln=3, 
            max_line_height=pdf_name.font_size)
      pdf_name.ln(line_height)

def combine_similar_columns(df, column_list):
    col_dict = {}
    for i in column_list:
        col_dict[i] = df[i].tolist()
    df_list = []
    for idx in range(len(df)):
        idx_list = []
        for col in col_dict.keys():
            if type(col_dict[col][idx]) == list:
                for item in col_dict[col][idx]:
                    idx_list.append(item)
            else:
                idx_list.append(col_dict[col][idx])
        idx_list = unique(remove_nan_from_list(idx_list))
        idx_list = ' '.join(idx_list)
        idx_list = idx_list.replace(",", ", ")
        df_list.append(idx_list)
    return df_list

def create_dt_col(df):
    dt_beta = df["dt_beta"]
    dt_omega = df["dt_omega"]
    dt_list = []
    for i in range(len(df)):
        if (dt_beta[i] or dt_omega[i]) == "positive":
            dt_list.append("+")
        else:
            dt_list.append("-")
    df["DT"] = dt_list
    return df

def create_dummy_data(df, col_name, fill_text):
    fill_list = []
    for i in range(len(df)):
        fill_list.append(fill_text)
    df[col_name] = fill_list
    return df

def join_pdfs(list, output_filename):
  output = PdfWriter()
  for pdf in list:
    reader = PdfReader(pdf)
    number_of_pages = len(reader.pages)

    pdf = PdfFileReader(open(pdf, "rb"))
    number_of_pages = len(pdf.pages)

    for i in range(number_of_pages):
      output.addPage(pdf.getPage(i))
  outputStream = open(output_filename, "wb")
  output.write(outputStream)
  outputStream.close()

def remove_nan_from_list(l):
    return [item for item in l if str(item) != 'nan']
  
def unique(list1): 
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

def _add_spaces(list, maxlen):
    head_list = []
    total_counts = []

    for string in list:
        string = str(string)
        if " " in string:
            new_str = string.split(" ")
            no_spaces = False
        else:
            new_str = string
            no_spaces = True
        new_line_count = []

        new_list = []
        for sub_str in new_str:
            new_sub = [sub_str[i:i+maxlen] for i in range(0, len(sub_str), maxlen)]
            new_sub = "\n".join(new_sub)
            new_list.append(new_sub)

        i = 0
        count = 0
        n_str = ' '.join(new_list)
        idx_list = []
        space_idx = []

        for i in range(len(n_str)):
            if n_str[i] == ' ':
                idx_list.append(i)
            elif n_str[i] == '\\'  and n_str[i+1]=='n':
                count = -1
                space_idx.append(i)
            else:
                count += 1
            if count == maxlen:
                if idx_list:
                    idx = _closest_value(idx_list, i)
                    if i-idx >=maxlen:
                        space_idx.append(i)
                        count = 0
                    else:
                        space_idx.append(idx)
                        count = i-idx
                else:
                    space_idx.append(i)
                    count = 0

        adding_list = []
        for space in range(len(space_idx)):
            adding_list.append(space)
        l = [x + y for x, y in zip(space_idx, adding_list)]
        new_list = " ".join(new_list)
        new_list = _convert(new_list)
        for i in l:
            new_list.insert(i, "\n")
        new_list = ''.join(new_list)
        if no_spaces:
            new_list = new_list.strip().replace(" ", "")
        head_list.append(new_list)
        new_line_count.append(new_list.count(' '))
        total_counts.append(max(new_line_count))
    return(head_list, new_line_count)

def _closest_value(input_list, input_value):
  arr = np.asarray(input_list)
  i = (np.abs(arr - input_value)).argmin()
  return arr[i]

def _convert(string):
    list1 = []
    list1[:0] = string
    return list1

def _flatten(l):
    return [item for sublist in l for item in sublist]

def _font_style(bold=False,underline=False,italic=False):
  font_style=""
  if bold:
    font_style+="B"
  if underline:
    font_style+="U"
  if italic:
    font_style+="I"
  return font_style

def _reformat(df, df_list=[],df_leng_list=[], max_line_len=40):
  df = df.copy()
  print(type(df))
  if not df_list:
      df_list = df.columns
  if not df_leng_list:
      df_leng_list=([max_line_len]*len(df))
  if len(df.columns) >=8:
    df = df.iloc[:, [0,1,2,3,4,5,6,7]]
    print("")
  max_height_list = []
  #Determine max number of rows of wrapped text per row of cells and create array of height multipliers
  print(df_list)
  for idx in range(len(df_list)):
      col = df_list[idx]
      col_list, max_vals = _add_spaces(df[col].tolist(), df_leng_list[idx])
      #print("max_val", max_vals, col_list)
      df[col] = col_list
      max_vals = max_vals[0]
      if max_vals == 0:
        max_vals = 1.1
      max_height_list.append(max_vals)
  
  df = df.replace('nan', '')
  df = df.replace(',',' ', regex=True)

  #Based on words per cell
  width_ratios = []
  #Adjust width ratio for cells with few words but many characters
  adjusting_values = []
  
  for col in df.columns:
      col_list = df[col].tolist()
      col_list.append(col)
      split_list = _flatten([row.split(" ") for row in _remove_nan_from_list(col_list)])
      max_word_len = max(split_list, key = len)
      width_ratios.append(len(max_word_len))
      adjusting_values.append(len(split_list))
  sum_ratio = sum(width_ratios)
  tenth = round(sum_ratio*.10)
  max_indices = list(np.array([index for index, item in enumerate(adjusting_values[2:]) if item == max(adjusting_values[2:])])+2)

  space = 0
  for i in range(len(width_ratios)):
      if i==0 or i==1:
          if width_ratios[i] - tenth >= 1:
              space+=(width_ratios[i] - tenth)
              width_ratios[i] = tenth
      elif i in max_indices:
          width_ratios[i]=width_ratios[i]+(space/len(max_indices))

  table = tuple(df.itertuples(index=False))
  return table, width_ratios, max_height_list

def _remove_nan_from_list(l):
    return [item for item in l if str(item) != 'nan']
  
def new_pdf(font="Times"):
  """Initiates new pdf file"""
  pdf = FPDF(orientation='P', unit='mm', format='letter')
  return pdf

__all__ = ["add_dendrogram_as_pdf", "add_page_header", \
"add_section_header", "add_paragraph", "add_table", "combine_similar_columns", \
"create_dt_col", "create_dummy_data", "join_pdfs", "remove_nan_from_list", "unique", \
new_pdf]# _add_spaces, _closest_value, _convert, _flatten, _font_style, _reformat, \
#_remove_nan_from_list]
