import shutil
from pathlib import Path

def create_output_directory(args)->None:
    directory = args.prefix
    path = Path(directory)
    print("Creating directory " + directory + " for output")
    path.mkdir(exist_ok=True)

    return path.as_posix()

def name_output_files(args, outputdir:str)->dict:
    files = {}
    files["testgenomebed"] = outputdir + "/genome." + args.assembly + ".bed"
    if not args.n_bedfile:
        files["testnbed"] = outputdir + "/nlocs." + args.assembly + ".bed"
    else:
        files["testnbed"] = None
    files["testnonnbed"] = outputdir + "/atgcseq." + args.assembly + ".bed"
    files["truthcovered"] = outputdir + "/" + args.assembly + ".benchcovered." + args.benchmark + ".bed"
    files["testmatcovered"] = outputdir + "/testmatcovered." + args.assembly + ".bed"
    files["testpatcovered"] = outputdir + "/testpatcovered." + args.assembly + ".bed"
    if args.variantfile is None:
        files["variantbed"] = outputdir + "/" + args.assembly + ".variantsinaligns." + args.benchmark + ".bed"
    else:
        files["variantbed"] = args.variantfile
    files["generalstatsfile"] = outputdir + "/" + args.assembly + ".generalstats.txt"
    files["mononucstatsfile"] = outputdir + "/" + args.assembly + ".mononucstats.txt"
    files["coveredmononucsfile"] = outputdir + "/" + args.assembly + ".coveredmononucs." + args.benchmark + ".bed"
    files["mononucswithvariantsfile"] = outputdir + "/" + args.assembly + ".mononucswithvariants." + args.benchmark + ".bed"
    files["bencherrortypebed"] = outputdir + "/" + args.assembly + ".errortype." + args.benchmark + ".bed"
    files["testerrortypebed"] = outputdir + "/errortype." + args.assembly + ".bed"
    files["coveredhetsitealleles"] = outputdir + "/" + args.assembly + ".coveredhetalleles." + args.benchmark + ".bed"

    return files
