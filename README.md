# CONSIDERACIONES GENERALES

Este programa tiene como función realizar la traducción a partir de una secuencia de nucleótidos, identificando los ORFs que se encuentran en la secuencia de DNA. Durante el proceso de traducción no se consideró el reverso complementario, por esta razón se sugiere utilizar una secuencia en dirección 5'- 3'.
Para identificar los ORFs, la secuencia de nucleotidos debe cumplir con las siguientes afirmaciones:

1. Iniciar con ATG.
2. Finalizar con las siguientes tripletas TGA, TAG, or TAA.
3. El número de bases totales del ORF debe ser múltiplo de 3.




# TRATAMIENTO DE EXCEPCIONES PRESENTES EN LA SECUENCIA DEL GEN 

La base de datos de proteinas fue construida para el genero Thysanoptera, y el programa fue diseñado con la secuencia fasta del gen mitrocondrial que codifica para la proteína cytochrome c oxidase subunit I, adicionalmente se calibró utilizando una secuencia modificada (nucleotide_prueba.txt) en la cual fueron incluidas algunas excepciones, entre ellas, minusculas (atgc), caracteres especiales (-----///****) y NNNNNN, para cada excepción se asumieron las siguientes acciones:

i) Si en la secuencia hay minúsculas estas serán reemplazadas por mayúsculas, en este caso, se utiliza el comando: 

   
   
   file = file.upper()

ii) Si en la secuencia de nucleótidos (file.fasta) hay caracteres diferentes a A T G C N, estos serán eliminados utilizando la siguiente linea de comandos, 



 

     
     for characters in file:
     
         if characters != "A" and characters != "T" and characters != "G" and characters != "C" and characters != "A" and characters != "N":
    
             file = file.replace(characters, "")
             
             
             
             
             
             
Con esta linea de comandos, se indica que si en la secuencia de nucleotidos  (file.fasta) hay caracteres diferentes a A, T, G, C, o N, estos deben reemplazarse por nada, es decir, seran eliminados, antes de realizar el proceso de traducción.  



iii) Si la secuencia de nucleotidos (file.fasta), incluye la letra N o secuencias de esta letra,  estas no serán eliminadas, y serán remplazadas por la letra ¨X¨ durante la traducción. Esta instrucción se da al definir la función de traducción. 


                        def translate(gen):
                        
                            protein_seq = ""
                            
                            for n in range(0, len(gen) +1, 3): # n = nucleotido
                            
                                codon = gen[n:n +3]
                                
                                if "N" in codon:
                                
                                    protein_seq += "X"
                                    
                                try:
                                
                                    protein_seq += codon_dictionary[codon]
                                    
                                except:
                                
                                    pass # Como en el diccionario no hay Ns, debe indicarse que ignore el error que se pueda generar.
                                    

                            return protein_seq
                            
                          
# IDENTIFICACIÓN DE ORFs

Para buscar los ORFs primero fue necesario identificar todos los codones de inicio presentes en la secuencia analizada, para esto se utilizaron los siguientes comandos


              """Encontrar sitios para iniciar la traducción."""
              
              # Generar arreglo con los sitios de inicio
              
              starts = []

              # Econtrar primer codon de inicio (remember, find() returns -1 if not found)
              
              i = sequence1.find('ATG')
              
              # Seguir buscando todos los codones de incio dentro de la secuencia

              while i >= 0:
              
                  starts.append(i)
                  
                  i = sequence1.find('ATG', i + 1)
                  
Con el objetivo de visualizar, que con esta linea de comandos se lograba encontrar todos los codones de inicio, se imprimieron los resultados con la 
siguiente instrucción

                  print(starts) # Este comando solo se utilizó para visualizar la ubicación de los codones de incio.
                  print(sequence1) # Este comando solo se utilizó para verificar que las posiciones generadas si correspondían al codon de inicio ATG
                  
               
<img width="1385" alt="image" src="https://user-images.githubusercontent.com/116923271/200409933-04c5137b-b76a-4659-92eb-7ed67ab2e546.png">


Después de verificar que la identifiación de los codones de incio funcionaba correctamente se procede define la función que se empleó para la busqueda de los codones de parada.


                    ## Definiendo función para encontrar el codon de parada.

                    print("Buscando codones de parada")

                    def find_stop(seq):

                        i = 0

                        while i < len(seq) - 2 and seq[i:i + 3] not in ('TAA', 'TAG', 'TGA'):
                            i += 3 # empiece a leer desde la posición cero hasta la última tripleta y sin incluir codones de parada TAA, TAG, TGA.

                     ## Si en encuentra el codon antes del final, imprima la posición del codon
                        if i < len(seq) -2:
                            stop = i + 3
                            ORF = seq[0:stop]

                            return ORF
                        else: # si no encuentra codon de parada pare lectura e imprima -1

                            return -1

Luego se verifica si la función para identificar los codones de parada corra correctamente,

      # # Verificando función para identificar codon de parada
      #
      # print("buscando codon de parada en la posicion:", find_stop(sequence1)) # solo se realizó para verificar si la función funcionaba correctamente
      # 
      # print(sequence1[0:find_stop(sequence1)]) # solo se realizó para verificar si la función para identificar el codon de parada funcionaba correctamente
      
# DEFINIENDO FUNCIÓN DE TRADUCCIÓN

Para realizar la traducción mediante python, primero fue necesario construir un diccionario para el código genético.

        # Creando diccionario con el código genético
                                            codon_dictionary = {"TTT": "F",
                                                                 "TTC": "F",
                                                                 "TTA": "L",
                                                                 "TTG": "L",
                                                                 "CTT": "L",
                                                                 "CTC": "L",
                                                                 "CTA": "L",
                                                                 "CTG": "L",
                                                                 "ATT": "I",
                                                                 "ATC": "I",
                                                                 "ATA": "I",
                                                                 "ATG": "M",
                                                                 "GTT": "V",
                                                                 "GTC": "V",
                                                                 "GTA": "V",
                                                                 "GTG": "V",
                                                                 "TCT": "S",
                                                                 "TCC": "S",
                                                                 "TCA": "S",
                                                                 "TCG": "S",
                                                                 "CCT": "P",
                                                                 "CCC": "P",
                                                                 "CCA": "P",
                                                                 "CCG": "P",
                                                                 "ACT": "T",
                                                                 "ACC": "T",
                                                                 "ACA": "T",
                                                                 "ACG": "T",
                                                                 "GCT": "A",
                                                                 "GCC": "A",
                                                                 "GCA": "A",
                                                                 "GCG": "A",
                                                                 "TAT": "Y",
                                                                 "TAC": "Y",
                                                                 "TAA": "*",
                                                                 "TAG": "*",
                                                                 "CAT": "H",
                                                                 "CAC": "H",
                                                                 "CAA": "Q",
                                                                 "CAG": "Q",
                                                                 "AAT": "N",
                                                                 "AAC": "N",
                                                                 "AAA": "K",
                                                                 "AAG": "K",
                                                                 "GAT": "D",
                                                                 "GAC": "D",
                                                                 "GAA": "E",
                                                                 "GAG": "E",
                                                                 "TGT": "C",
                                                                 "TGC": "C",
                                                                 "TGA": "*",
                                                                 "TGG": "W",
                                                                 "CGT": "R",
                                                                 "CGC": "R",
                                                                 "CGA": "R",
                                                                 "CGG": "R",
                                                                 "AGT": "S",
                                                                 "AGC": "S",
                                                                 "AGA": "R",
                                                                 "AGG": "R",
                                                                 "GGT": "G",
                                                                 "GGC": "G",
                                                                 "GGA": "G",
                                                                 "GGG": "G"}


La función traducción, se define con la siguiente linea de comandos.

              def translate(gen):
                  protein_seq = ""
                  for n in range(0, len(gen) +1, 3): # n = nucleotido
                      codon = gen[n:n +3]
                      if "N" in codon:
                          protein_seq += "X"
                      try:
                          protein_seq += codon_dictionary[codon]
                      except:
                          pass #como en el diccionario no hay Ns, debe indicarse que ignore el error que se genere.

                  return protein_seq

Para poder realizar este paso de la traducción, fue necesario incluir la convención para los codones que contengan ´Ns´ mencionada al inicio de este documento, en la cual se indica que en caso de presentarse ´Ns´ dentro o en las tripletas estas seran remplazadas por la letra "X", que representa cualquier aminoácido.

Con las funciones definidas, se procede a realizar el proceso de traducción.

1. Generación de preORFs.

Para generar el listado de todos los preORFs es necesario unir las secuencias generadas para cada codon de inicio, para lograrlo se corre el siguiente comando

      listpreORF = []

      for scodon in starts:
          preORF = sequence1[scodon:]
          listpreORF.append(preORF) # Agregar todos los preORFs identificados a la lista

      print(listpreORF) # Este comando solo se utilizó para verificar que todos los preORFs se agregaron a la lista.
,
 
<img width="946" alt="image" src="https://user-images.githubusercontent.com/116923271/200424199-61561afe-e8dd-4de5-a2fe-5840585c6f7c.png">

Después de verificar que todos los preORFs se encontraban en una lista, se procede con la traducción y se genera el archivo fasta para las proteinas encontradas.

2. Traducción y generación de archivo fasta de las proteinas encontradas.

La traducción se realizó con la siguiente línea de comandos, y llamando la función (translate) para la traducción y la función para encontrar el codón de parada (find_stop) definidas previamente.

                  f = open("proteinas.fasta", "a") # Archivo para escribir las proteinas encontradas

                  s = 0  # Asignando contador para la iniciar la lectura de la secuencia

                  # Generando formato .fasta

                  print("Realizando traducción")

                  for line in listpreORF:
                  
                      s += 1 #contando
                      
                      protein_seq = translate(find_stop(line)) #llamando funciones translate y find_stop
                      
                      f.write(f'>{s}\n') #enumerando proteinas
                      
                      f.write(f'{protein_seq}\n')

                      print(protein_seq) # imprimiendo proteinas
                      
                  f.close()

El resultado arrojado por esta línea de comandos es un archivo llamado ¨proteinas.fasta¨, el cual será utilizado para correr el blastp.

<img width="492" alt="image" src="https://user-images.githubusercontent.com/116923271/200433954-4c056562-55bb-4261-ab51-022a4e4beee2.png">

3. Identificación de proteinas

En este punto ya se procede a correr el archivo fasta de las proteinas encontradas frente a una base de datos de proteinas para Thrips. Pero antes se debe crear la base de datos con secuencias de referencia descargadas del NCBI (https://www.ncbi.nlm.nih.gov/protein/?term=thysanoptera), utilizando la siguiente linea de comandos.

                Print('Creando Base de Datos proteins_Thripsdb.fasta')

                os.system('makeblastdb -dbtype prot -in proteins_Thripsdb.fasta')  # Con este comando se puede crear desde la terminal la base de datos usando el comando de la rutina

Finalmente, se corre el blastp de las proteinas encontradas frente a la base de datos del género thysanoptera, el comando utilizado es el siguiente

                comando_blast = f'blastp -db proteins_Thripsdb.fasta -query proteinas.fasta -outfmt "6 qseqid stitle pident evalue qcovs" -out resultados.tsv'

                os.system(comando_blast) # con este comando se puede correr desde la terminal el blastp usando el comando de la rutina

 
 Los resultados se guaradaran en el archivo ¨resultados.tsv¨, donde se podrán observar los nombres de las proteinas encontradas en el archivo analizado con el programa.
 
 <img width="872" alt="image" src="https://user-images.githubusercontent.com/116923271/200439596-d7368503-4f94-437b-85e8-f16e4cda9c53.png">
 
 
 #ANEXOS
 
 Archivos utilizados en la calibración del programa
 
 1. Archivo de modificado con excepciones
 [nucleotide_prueba.txt](https://github.com/mmhiguit/Translation/files/9956379/nucleotide_prueba.txt)
 2. Secuencia original sin modificar
 [nucleotide_Thrips.txt](https://github.com/mmhiguit/Translation/files/9956387/nucleotide_Thrips.txt)
 3. Secuencias de referencia del NCBI de proteinas utilizadas para la creación de la base de datos.
 [Copy_proteins_Thripsdb.txt.zip](https://github.com/mmhiguit/Translation/files/9956409/Copy_proteins_Thripsdb.txt.zip)


 

 


















