# ClonalTree

**Reconstructing B cell lineage trees with minimum spanning tree and genotype abundances**

ClonalTree is a new algorithm to reconstruct BCR lineage trees that incorporates genotype abundance into a minimum spanning tree to infer maximum parsimony trees.

**CONTACT**  
  E-mail: 
  juliana.silva_bernardes@sorbonne-universite.fr 
  
## Inputs
 
  * The alignment of clonally related BCR heavy chain sequences in FASTA format. The naive B cell heavy chain sequence must be included in the alignment, named as "naive". There are ways to provide the genotype abundance for each sequence :
    * By the sequence ID repetition in the alignment, or
    * By integrating the abundance of each sequence in its ID, for instance, the sequence Seq1 with an abundance of 6 will have the following ID:
      >Seq1@6
  * See [example input files](https://github.com/julibinho/ClonalTree/tree/main/Examples/input)

## Outputs

  * ClonalTree returns:

    - [repertoire_name].nk : the reconstructed BCR lineage tree in [newick format](https://en.wikipedia.org/wiki/Newick_format) 

    - [repertoire_name].nk.csv :  a table in csv format, containing the parent relationship and cost.
  * See [example output files](https://github.com/julibinho/ClonalTree/tree/main/Examples/output)

     
      
## Requirements 

  * numpy :
      ```
      $ conda install numpy
      ```
      or 
      ```
      $ pip install numpy
      ```

  * Biopython
      ```
      $ pip install biopython
      ```

  * ete3 :
      ```
      $ pip install ete3
      ```

  * networkx (for the evaluation) :
      ```
      $ pip install networkx
      ```


## Using ClonalTree 
   The command line for launching the ClonalTree is:

  ```
  $ python clonalTree.py  -i [seq_alignment_file] -o [output_newick_file] [...options]

  ```
### required arguments 
  * [seq_alignment_file] is alignment of clonally related BCR heavy chain sequences in FASTA format,
  * [output_newick_file] is output file name

### optional arguments [...options]

  * -a 1, if considering abundance, otherwise -a 0
  * -r 1, if performing revision, otherwise -r 0
  * -t 1, if performing trimming tree, otherwise -t 0


  For instance the following command can be run in the src/ folder:
  ```
  $ python src/clonalTree.py  -i Examples/input/simulation200.fasta -o Examples/output/clonalTree.abRT.nk -a 1 -r 1 -t 1
  ```
                      
  Output files will be placed as such:
  ```
  ~Examples/output/[seq_alignment_file].nk
                  [seq_alignment_file].nk.csv
 ```
 [seq_alignment_file] is the multiple sequence alignement of clonaly-related B cell receptor heavly chain seuqneces.

## License, Patches, and Ongoing Developements

  * The program is distributed under the CeCILL licence 
  * [Feature requests and open issues](https://github.com/julibinho/ClonalTree/issues).
 
 
## Reference
Nika Abdollahi, Anne Langlois De Septenville,  Frederic Davi and Juliana S. Bernardes. Reconstructing B cell lineage trees with minimum spanning tree and genotype abundances. BMC Bioinformatics (2023): 24(1) 70
