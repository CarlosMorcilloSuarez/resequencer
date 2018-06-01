        resequencer.py
        
            Generates a simulated fastq file with reads from a NGS experiment
            from a reference genome in fasta format
            
            Use
                resequencer.py [-r,--reference] genome.fa [options]
                
                
        fastResequencer.py   

            It's a faster and less memory consuming version of resequencer.py but
            doesn't allow for SVs (structural variants) introduction.
            
            Use
                fastResequencer.py [-r,--reference] genome.fa [options]
