# Attention, ce code s'adresse principalement à des fichiers dont la structure ressemble à ceci (il y a un # en trop sur chaque ligne, c'est nécessaire pour tout mettre en commentaires): 
# # STOCKHOLM 1.0
# #=GF ID   DUF1798
# #=GF AC   PF08807.15
# #=GF DE   Bacterial domain of unknown function (DUF1798)
# #=GF GA   25; 25;
# #=GF TP   Domain
# #=GF ML   108
# //
# # STOCKHOLM 1.0
# #=GF ID   DUF3642
# #=GF AC   PF12182.13
# #=GF DE   Bacterial lipoprotein
# #=GF GA   22.7; 22.7;
# #=GF TP   Domain
# #=GF ML   93
# #=GF CL   CL0116
# //
# # STOCKHOLM 1.0
# #=GF ID   DUF4128
# #=GF AC   PF13554.11
# #=GF DE   Bacteriophage related domain of unknown function
# #=GF GA   27; 27;
# #=GF TP   Family
# #=GF ML   127
# #=GF CL   CL0691
# //

# Ce qui est considéré comme étant un "motif" dans le programme est un ensemble de ligne qui se répète du type : 

# # STOCKHOLM 1.0
# #=GF ID   exemple
# #=GF AC   PF99999
# #=GF DE   Exemple exemple exemple...
# #=GF GA   00; 00;
# #=GF TP   Exemple
# #=GF ML   999
# #=GF CL   EX0000
# //

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

def main():

    parser = argparse.ArgumentParser(description='read a dat file with a specific structure and filter in order to have only the elements of interest in a new file created for this purpose')
    parser.add_argument('dat_file', type=str, help='dat file of interest which you want to filter')
    parser.add_argument('out_file',type=str, help='output file that will be created')
    parser.add_argument('keyword', type=str, help='the word that will serve as filter, watch out to lowercases and uppercases')

    args = parser.parse_args()
    dat_file = args.dat_file
    out_file = args.out_file
    keyword = args.keyword

    motif_extractor(dat_file, out_file, keyword)

if __name__ == "__main__":
    main()
