import os
import glob
import re
import importlib.resources

def plot_benchmark_align_coverage(assemblyname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('BenchCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " + resourcedir
        print(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_testassembly_align_coverage(assemblyname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('TestCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + outputdir + " " + resourcedir
        print(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_mononuc_accuracy(assemblyname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('MononucAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + outputdir + " " + resourcedir
        print(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_svcluster_align_plots(assemblyname:str, benchname:str, outputdir:str, resourcedir:str, refobj):
    rfile_res = importlib.resources.files("q100bench").joinpath('PlotChromAligns.R')

    chromalignbedfiles = glob.glob(outputdir + "/*.clusters.bed")
    returnvalues = []
    with importlib.resources.as_file(rfile_res) as rfile:
        for chrombed in chromalignbedfiles:
            chromosome = chrombed.replace(".clusters.bed", "")
            chromosome = re.sub(r".*/*clustered_aligns\.", "", chromosome)
            chromlength = refobj.get_reference_length(chromosome)
            plotcommand = "Rscript " + str(rfile) + " " + chrombed + " " + assemblyname + " " + benchname + " " + outputdir + " " + resourcedir + " " + str(chromlength)
            print(plotcommand)
            returnvalues.append(os.system(plotcommand))
    return returnvalues

