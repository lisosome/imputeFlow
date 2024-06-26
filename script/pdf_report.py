#!/usr/bin/env python3

import pandas as pd
from fpdf import FPDF
import argparse
from glob import glob
import sys
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mode", help="Select mode of execution", choices=["CHUNK", "CHR"])
    parser.add_argument("-r", "--chr", help="Current chromosome")
    parser.add_argument("-b", "--base", help = "Stats basefolder")
    parser.add_argument("-o", "--outfile", help="Name of the output pdf report")
    args=parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args


def get_chunk_n(stats_folder):
    """get_chunk_n

    Args:
        stats_folder (str): Path in which chinks report are stored
    The function will search for a type of report generated per chunk in order to determine the total chunk number per chromosome
    """
    chrom=os.path.basename(os.path.dirname(stats_folder))
    filename=f"""{chrom}_*_impute_summary.png"""
    file_to_search=os.path.join(stats_folder, filename)
    lst = glob(file_to_search)
    tot_chunk_n = len(lst)
    return tot_chunk_n


def pdf_report_chunk(current_chr,stat_base_folder,chunk_n,outfile):

    class PDF(FPDF):
        def footer(self):
            # Position at 1.5 cm from bottom
            self.set_y(-15)
            # helvetica italic 8
            self.set_font('Times', 'I', 8)
            # Page number
            self.cell(0, 10, 'Page ' + str(self.page_no()) + '/{nb}', 0, 0, 'C')


    #we need to build the template for the single chunk, than collect all chunks data together
    # define the current page elements
    # current_chr='2' # this is something we need to pass as arguments of our function
    # stat_base_folder="/home/cocca/analyses/imputation/20210613/MOLISANI/07.stats/2/CHUNKS" #this is another argument needed to get the rest of the data

    #open the PDF documents and start the page numbering
    pdf = PDF()
    # pdf.alias_nb_pages()

    for chunk in list(range(1,int(chunk_n)+1)):
        #current_chunk="{:02d}".format(chunk)
        current_chunk=str(chunk)
        #define files needed for the pdf generation
        current_info_af=stat_base_folder+'/'+current_chr+'_'+ current_chunk +'_impute_summary.png'
        current_manhattan=stat_base_folder+'/'+current_chr+'_'+ current_chunk +'_impute_manhattan.png'
        current_chunk_stats_by_maf=stat_base_folder+'/'+current_chr+'_'+ current_chunk +'_impute_summary_by_maf.csv'
        current_chunk_stats_by_maf_by_info=stat_base_folder+'/'+current_chr+'_'+ current_chunk +'_impute_summary_by_maf_by_info.csv'
        #create the new page for the current chunk
        pdf.add_page()
        pdf.set_font('Times', 'B', 16)
        #add plots
        page_header="Info score and allele frequency report for chr" + current_chr + " chunk " + current_chunk
        pdf.cell(0, 10, page_header, ln=1, align='C')
        pdf.image(current_info_af,w=200)
        pdf.ln(1)
        pdf.image(current_manhattan,w=190)
        #now add table with report for each chunk
        #this is something generated by pandas, so should be easy to read in
        df_stats_by_maf=pd.read_csv(current_chunk_stats_by_maf,index_col=0)
        df_stats_by_maf_by_info=pd.read_csv(current_chunk_stats_by_maf_by_info,names=["INFO_CLASS","MAF_CLASS","count","mean","std","min","25%","50%","75%","max"], skiprows=3)

        stats_by_maf = list(df_stats_by_maf.to_records(index=True))
        stats_by_maf_header=list(df_stats_by_maf.columns)
        stats_by_maf_header.insert(0,'')
        stats_by_maf.insert(0,tuple(stats_by_maf_header))
        stats_by_maf_by_info = list(df_stats_by_maf_by_info.to_records(index=False))
        stats_by_maf_by_info.insert(0,tuple(df_stats_by_maf_by_info.columns))

        stats_by_maf_header='Stats of INFO scores by MAF classes'
        stats_by_maf_by_info_header='Stats of INFO scores by MAF classes and by INFO score classes'
        #add page with tables
        pdf.add_page()
        pdf.set_font("Times", size=10)
        pdf.cell(0, 10, stats_by_maf_header, ln=1, align='C')
        line_height = pdf.font_size * 2.5
        col_width = pdf.epw / 10  # distribute content evenly
        for row in stats_by_maf:
            for datum in row:
                if isinstance(datum,str) :
                    pdf.multi_cell(col_width, line_height, datum, border=1, ln=3, max_line_height=pdf.font_size,align='C')
                else:
                    pdf.multi_cell(col_width, line_height, str(round(datum,5)), border=1, ln=3, max_line_height=pdf.font_size,align='C')
            pdf.ln(line_height)

        pdf.ln(1)
        pdf.cell(0, 10, stats_by_maf_by_info_header, ln=1, align='C')
        for row in stats_by_maf_by_info:
            for datum in row:
                if isinstance(datum,str) :
                    pdf.multi_cell(col_width, line_height, datum, border=1, ln=3, max_line_height=pdf.font_size,align='C')
                else:
                    pdf.multi_cell(col_width, line_height, str(round(datum,5)), border=1, ln=3, max_line_height=pdf.font_size,align='C')
            pdf.ln(line_height)

    pdf.output(outfile)


def pdf_report_chr(current_chr,stat_base_folder,outfile):

    class PDF(FPDF):
        def footer(self):
            # Position at 1.5 cm from bottom
            self.set_y(-15)
            # helvetica italic 8
            self.set_font('Times', 'I', 8)
            # Page number
            self.cell(0, 10, 'Page ' + str(self.page_no()) + '/{nb}', 0, 0, 'C')


    #we need to build the template for the single chunk, than collect all chunks data together
    # define the current page elements
    # current_chr='2' # this is something we need to pass as arguments of our function
    # stat_base_folder="/home/cocca/analyses/imputation/20210613/MOLISANI/07.stats/2/CHUNKS" #this is another argument needed to get the rest of the data

    #open the PDF documents and start the page numbering
    pdf = PDF()

    current_info_af=stat_base_folder+'/'+current_chr+'_impute_summary.png'
    current_manhattan=stat_base_folder+'/'+current_chr+'_impute_manhattan.png'
    current_chunk_stats_by_maf=stat_base_folder+'/'+current_chr+'_impute_summary_by_maf.csv'
    current_chunk_stats_by_maf_by_info=stat_base_folder+'/'+current_chr+'_impute_summary_by_maf_by_info.csv'
    #create the new page for the current chunk
    pdf.add_page()
    pdf.set_font('Times', 'B', 16)
    #add plots
    page_header="Info score and allele frequency report for chr" + current_chr
    pdf.cell(0, 10, page_header, ln=1, align='C')
    pdf.image(current_info_af,w=200)
    pdf.ln(1)
    pdf.image(current_manhattan,w=190)
    #now add table with report for each chunk
    #this is something generated by pandas, so should be easy to read in
    df_stats_by_maf=pd.read_csv(current_chunk_stats_by_maf,index_col=0)
    df_stats_by_maf_by_info=pd.read_csv(current_chunk_stats_by_maf_by_info,names=["INFO_CLASS","MAF_CLASS","count","mean","std","min","25%","50%","75%","max"], skiprows=3)

    stats_by_maf = list(df_stats_by_maf.to_records(index=True))
    stats_by_maf_header=list(df_stats_by_maf.columns)
    stats_by_maf_header.insert(0,'')
    stats_by_maf.insert(0,tuple(stats_by_maf_header))
    stats_by_maf_by_info = list(df_stats_by_maf_by_info.to_records(index=False))
    stats_by_maf_by_info.insert(0,tuple(df_stats_by_maf_by_info.columns))

    stats_by_maf_header='Stats of INFO scores by MAF classes'
    stats_by_maf_by_info_header='Stats of INFO scores by MAF classes and by INFO score classes'
    #add page with tables
    pdf.add_page()
    pdf.set_font("Times", size=10)
    pdf.cell(0, 10, stats_by_maf_header, ln=1, align='C')
    line_height = pdf.font_size * 2.5
    col_width = pdf.epw / 10  # distribute content evenly
    for row in stats_by_maf:
        for datum in row:
            if isinstance(datum,str) :
                pdf.multi_cell(col_width, line_height, datum, border=1, ln=3, max_line_height=pdf.font_size,align='C')
            else:
                pdf.multi_cell(col_width, line_height, str(round(datum,5)), border=1, ln=3, max_line_height=pdf.font_size,align='C')
        pdf.ln(line_height)

    pdf.ln(1)
    pdf.cell(0, 10, stats_by_maf_by_info_header, ln=1, align='C')
    for row in stats_by_maf_by_info:
        for datum in row:
            if isinstance(datum,str) :
                pdf.multi_cell(col_width, line_height, datum, border=1, ln=3, max_line_height=pdf.font_size,align='C')
            else:
                pdf.multi_cell(col_width, line_height, str(round(datum,5)), border=1, ln=3, max_line_height=pdf.font_size,align='C')
        pdf.ln(line_height)

    pdf.output(outfile)



def main():
    
    args = get_args()
    mode = args.mode
    chrom = args.chr
    stat_folder = args.base
    outfile = args.outfile

    if mode == "CHUNK":
        chunk_n = get_chunk_n(stat_folder)
        pdf_report_chunk(chrom,stat_folder,chunk_n,outfile)
    else:
        pdf_report_chr(chrom,stat_folder,outfile)


if __name__ == '__main__':
    main()
