import os
import importlib.resources

def plot_benchmark_align_coverage(assemblyname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('BenchCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + outputdir + " " + resourcedir
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_testassembly_align_coverage(assemblyname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('TestCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + outputdir + " " + resourcedir
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_mononuc_accuracy(assemblyname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('MononucAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + outputdir + " " + resourcedir
        returnvalue = os.system(plotcommand)
    return returnvalue

