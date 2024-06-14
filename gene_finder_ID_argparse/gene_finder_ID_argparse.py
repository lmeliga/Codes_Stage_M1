import re # module indispensable pour les séparateurs plus bas dans le code 
import argparse
                 
def gene_finder_ID(gff_file, gene_wanted): # la fonction prend en arguments le fichier gff avec les annotations (le GCA) et le gène recherché (attention à l'orthographe)
    count_genes = 0 # création et initialisation d'un compteur pour compter le nombre de gènes recherchés présents dans le gff au total
    gene = [] # création d'une liste vide qui contiendra les lignes des gènes visés
    with open(gff_file, 'r') as gff : # on ouvre le fichier gff en question (on le renomme "gff" juste dans la fonction pour faciliter l'écriture)
         for line in gff: # Pour chaque ligne du fichier, on effectue la commande en dessous
            if not line.startswith("#") : # Si la ligne ne commence PAS par un #...
                colonne = line.strip().split("\t") # On définit une variable colonne. On divise la ligne en colonnes en utilisant la tabulation en guise de séparateur et on retire tout potentiel espace en début et en fin de ligne 
                if len(colonne) < 9 or colonne[2] != 'CDS': # Si la ligne possède moins de 9 colonnes ou qu'il n'y a pas écrit "CDS" à la troisième colonne (normalement il y a écrit "gene" ou "CDS")
                    continue # Si une des deux conditions précédentes n'est pas respectée, on passe à la ligne suivante
                qualifiers = "\t".join(colonne[8:]) # On fusionne toutes les colonnes à partir de la 9ème jusqu'à la dernière (peu importe le nombre) en une seule chaine de caractères en utilisant la tabulation comme séparateur 
                second_part = re.split(r'[ =\t:;()]+', qualifiers) # définit la "seconde partie" de la ligne, et divise uniquement cette dernière (grâce à qualifiers) selon les séparateurs qui sont l'espace, =, la tabulation, :, ;, (), et les nombres de 0 à 9. Le "+" indique qu'il peut y avoir plusieurs de chacun de ces séparateurs
                if gene_wanted in second_part: # Si le gène recherché est dans la "seconde partie"...
                    gene.append(line.strip()) # ...on ajoute la ligne dans la liste nommée "gene"
                    count_genes += 1 # ... et on incrémente le compteur de 1
    return gene, count_genes

def main() :

    parser = argparse.ArgumentParser(description='Use a GFF file with annotations to find the wanted gene from his ID')
    parser.add_argument('gene_wanted', type=str, help='Enter the ID of the gene that you are looking for. It should look like this : L_N. L is a one or several letters and N is one or several numbers')
    parser.add_argument('gff_file', type=str,help='gff file with the annotations')
    args = parser.parse_args()
    gene_wanted = args.gene_wanted
    gff_file = args.gff_file
    genes, number = gene_finder_ID(gff_file, gene_wanted)

    upgraded_genes = "\n\n".join(genes) # On prend la liste "genes" ("genes" et non "gene") et on crée une chaine de caractères avec "join" et entre chaque élement, on intègre un saut de ligne

    print("Gènes trouvés :")
    print(upgraded_genes)
    print("\nNombre de gènes trouvés :", number)

if __name__ == "__main__":
    main()