#### IMPORTS ####
import os
import textwrap
## Written by Roshan, Bereket, Paul

#### SEQUENCES ####
seqD = {
    "TBXT1": """MSSPGTESAGKSLQYRVDHLLSAVENELQAGSEKGDPTERELRVGLEESELWLRFKELTNEMIVTKNGRRMFPVLKVNVSGLDPNAMYSFLLDFVAADNHRWKYVNGEWVPGGKPEPQAPSCVYIHPDSPNFGAHWMKAPVSFSKVKLTNKLNGGGQIMLNSLHKYEPRIHIVRVGGPQRMITSHCFPETQFIAVTAYQNEEITALKIKYNPFAKAFLDAKERSDHKEMMEEPGDSQQPGYSQWGWLLPGTSTLCPPANPHPQFGGALSLPSTHSCDRYPTLRSHRSSPYPSPYAHRNNSPTYSDNSPACLSMLQSHDNWSSLGMPAHPSMLPVSHNASPPTSSSQYPSLWSVSNGAVTPGSQAAAVSNGLGAQFFRGSPAHYTPLTHPVSAPSSSGSPLYEGAAAATDIVDSQYDAAAQGRLIASWTPVSPPSM""",
    "TBXTG177D": """MSSPGTESAGKSLQYRVDHLLSAVENELQAGSEKGDPTERELRVGLEESELWLRFKELTNEMIVTKNGRRMFPVLKVNVSGLDPNAMYSFLLDFVAADNHRWKYVNGEWVPGGKPEPQAPSCVYIHPDSPNFGAHWMKAPVSFSKVKLTNKLNGGGQIMLNSLHKYEPRIHIVRVGDPQRMITSHCFPETQFIAVTAYQNEEITALKIKYNPFAKAFLDAKERSDHKEMMEEPGDSQQPGYSQWGWLLPGTSTLCPPANPHPQFGGALSLPSTHSCDRYPTLRSHRSSPYPSPYAHRNNSPTYSDNSPACLSMLQSHDNWSSLGMPAHPSMLPVSHNASPPTSSSQYPSLWSVSNGAVTPGSQAAAVSNGLGAQFFRGSPAHYTPLTHPVSAPSSSGSPLYEGAAAATDIVDSQYDAAAQGRLIASWTPVSPPSM"""
}

#### COMPLEX DEFINITIONS ####
complexes = {
    "complex1": ["TBXT1"],
    "complex2": ["TBXTG177D"]
}

#### OUTPUT DIR ####

def createDir(dirName = "tbxt"):
    wd = os.getcwd()
    outputDir = os.path.join(wd, dirName)
    print(f"Creating yaml directory {outputDir}")
    os.makedirs(outputDir, exist_ok=True)
    return outputDir

#### YAML GENERATOR ####
def make_yaml(complex_name, chains, outputDir):
    """Write a YAML file for a given complex with its protein chain(s)."""
    lines = []
    lines.append("sequences:")

    # Add each protein chain
    for i, protein_id in enumerate(chains):
        seq = seqD[protein_id]
        chain_id = chr(65 + i)
        lines.append(f"  - protein:")
        lines.append(f"      id: {chain_id}")
        lines.append(f"      sequence: |")
        for chunk in textwrap.wrap(seq, width=80):
            lines.append(f"        {chunk}")
        lines.append(f"      cyclic: false")

    # Define complex composition
    lines.append("complexes:")
    protein_ids = [chr(65 + i) for i in range(len(chains))]
    lines.append(f"  - proteins: [{', '.join(protein_ids)}]")

    # Write YAML file
    filename = f"{complex_name}.yaml"
    filepath = os.path.join(outputDir, filename)
    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"✅ Wrote {filepath}")

#### MAIN ####
if __name__ == "__main__":
    outputDir = createDir()
    for complex_name, chains in complexes.items():
        make_yaml(complex_name, chains, outputDir = outputDir)