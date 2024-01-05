import os

def plot_benchmark_align_coverage(assemblyname:str, outputdir:str, resourcedir:str):
    plotcommand = "Rscript /data/Phillippy/tools/q100bench/q100bench/BenchCoveragePlot.R " + assemblyname + " " + outputdir + " " + resourcedir
    returnvalue = os.system(plotcommand)
    return returnvalue

def plot_testassembly_align_coverage(assemblyname:str, outputdir:str, resourcedir:str):
    plotcommand = "Rscript /data/Phillippy/projects/HG002_diploid/benchmarking/software/TestCoveragePlot.R " + assemblyname + " " + outputdir + " " + resourcedir
    returnvalue = os.system(plotcommand)
    return returnvalue

def plot_mononuc_accuracy(assemblyname:str, outputdir:str, resourcedir:str):
    plotcommand = "Rscript /data/Phillippy/projects/HG002_diploid/benchmarking/software/MononucAccuracy.R " + assemblyname + " " + outputdir + " " + resourcedir
    returnvalue = os.system(plotcommand)
    return returnvalue

