import os
import glob
import re
import logging
import importlib.resources

logger = logging.getLogger(__name__)

def plot_benchmark_align_coverage(assemblyname:str, benchname:str, outputdir:str, benchparams:dict):
    rfile_res = importlib.resources.files("q100bench").joinpath('BenchCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        genomefile = benchparams["genomeregions"]
        nlocfile = benchparams["nstretchregions"]
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " + genomefile + " " + nlocfile
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_testassembly_align_coverage(assemblyname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('TestCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + outputdir + " " + resourcedir + " " + benchname
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_mononuc_accuracy(assemblyname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('MononucAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " + resourcedir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_qv_score_concordance(assemblyname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('PlotAssemblyQualValueAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " + resourcedir
        logger.info(plotcommand)
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
            logger.debug(plotcommand)
            returnvalues.append(os.system(plotcommand))
    return returnvalues

def plot_assembly_error_stats(assemblyname:str, genomename:str, outputdir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('IndelLengthPlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + genomename + " " + outputdir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_read_error_stats(readsetname:str, genomename:str, outputdir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('IndelLengthPlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + readsetname + " " + genomename + " " + outputdir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_read_mononuc_stats(readsetname:str, genomename:str, outputdir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('ReadMononucAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + readsetname + " " + genomename + " " + outputdir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_read_qv_score_concordance(readsetname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("q100bench").joinpath('PlotAssemblyQualValueAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + readsetname + " " + benchname + " " + outputdir + " " + resourcedir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

