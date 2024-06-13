import re
import argparse

def motif_extractor(input_file,  output_file, keyword):
    with open(input_file, 'r') as infile, open(output_file,'w') as outfile: # On ouvre le fichier d'entrée en mode lecture et le fichier de sortie en mode écriture (s'il n'existe pas dans le répertoire, il est créé). Avec "with", on s'assure que les fichiers sont correctement fermés après utilisation.
        motif = [] # Liste vide pour stocker les lignes du motif en cours de traitement 
        inside_motif = False # Booléen qui permet de dire si le programme est actuellement à l'intérieur d'un motif ou non. Par logique, ici, on l'initialise à faux. 

        for line in infile: # on lit le fichier d'entrée ligne par ligne 
            if line.startswith('# STOCKHOLM 1.0'): # Vérifie que la ligne actuelle commence par ça : "# STOCKHOLM 1.0". Dans ce cas, on est au début d'un nouveau motif.
                if motif: # Vérifie si la liste contient déjà des lignes (ce qui signifie qu'un motif présent précédent est toujours en cours de traitement par le programme).
                    if any(keyword in item for item in motif): # Si le motif en cours contient des lignes, vérifie si l'une de ces lignes contient le mot-clé (keyword). Si c'est le cas, écrit le motif dans le fichier de sortie. POur faire simple, cette ligne signifie : "Si le mot-clé est présent dans au moins une des lignes du motif..." 
                        outfile.write(''.join(motif))
                motif = [line] # nouveau motif initialisé en ajoutant la ligne actuelle à la liste de 'motif'
                inside_motif = True # nous sommes dans un nouveau motif donc true 
            elif line.startswith('//'): # Vérifie si la ligne commence par "//", ce qui indique la fin d'un motif
                motif.append(line)
                if any(keyword in item for item in motif): # Vérifie si l'une des lignes du motif en cours contient le mot-clé (keyword). Si c'est le cas, écrit le motif dans le fichier de sortie.
                    outfile.write(''.join(motif))
                motif = [] # Réinitialise la liste et définit inside_motif à false pour indiquer que nous ne sommes plus dans un motif 
                inside_motif = False
            elif inside_motif : # Vérifie que nous sommes au sein d'un motif 
                motif.append(line) # Ajoute la ligne actuelle au motif en cours 

        if motif and any(keyword in item for item in motif): # Après la fin de la boucle (c'est-à-dire après avoir lu toutes les lignes du fichier), vérifie si la liste motif contient des lignes (au cas où le fichier ne se termine pas par //). Si le dernier motif contient le mot-clé (keyword), l'écrit dans le fichier de sortie.
            outfile.write(''.join(motif)) 

def database_limiter(input_file, acc_list, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        motif = []
        inside_motif = False 

        for line in infile:
            if line.startswith('HMMER3/f'):
                if motif:
                    if any(any(acc in item for item in motif) for acc in acc_list):
                        outfile.write(''.join(motif))
                motif = [line]
                inside_motif = True
            elif line.startswith('//'):
                motif.append(line)
                if any(any(acc in item for item in motif) for acc in acc_list):
                    outfile.write(''.join(motif))
                motif = []
                inside_motif = False
            elif inside_motif:
                motif.append(line)

        if motif and any(any(acc in item for item in motif) for acc in acc_list):
            outfile.write(''.join(motif))

def main():    

    parser = argparse.ArgumentParser(description='read a file (containing the hmm patterns) with a specific structure and filter in order to have only the elements of interest in a new file created for this purpose')
    parser.add_argument('dat_file', type=str,help='dat file of interest which you want to filter')
    parser.add_argument('out_file',type=str, help='first output file that will be created')
    parser.add_argument('dat_file2', type=str,help='file of interest containing the HMM profiles')
    parser.add_argument('out_file2',type=str, help='last output file that will be created')
    parser.add_argument('keyword', type=str, help='the word that will serve as filter, watch out to lowercases and uppercases')

    args = parser.parse_args()
    dat_file = args.dat_file
    out_file = args.out_file
    dat_file2 = args.dat_file2
    out_file2 = args.out_file2
    keyword = args.keyword

    acc_list = []
    motif_extractor(dat_file, out_file, keyword)

    with open(out_file, 'r') as file:
        for line in file:
            if line.startswith("# STOCKHOLM 1.0"):
                continue
            if line.startswith("#=GF AC"):
                colonne = re.split(r'[ \t\n]+', line)
                acc_list.append(colonne[2])

    database_limiter(dat_file2, acc_list, out_file2)

if __name__ == "__main__":
    main()
