# rsFQC
FastQ Quality Control in rust - rustFastQC with in-terminal graphs thanks to [textplots-rs](https://github.com/loony-bean/textplots-rs/)  
[![Rust](https://github.com/a-xavier/rsFQC/actions/workflows/rust.yml/badge.svg)](https://github.com/a-xavier/rsFQC/actions/workflows/rust.yml)

### Currently testing:  
- Size distribution
- Quality distributions
- Duplications

### Instructions  

[Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) and rust needs to be installed first.  

```git clone https://github.com/a-xavier/rsFQC```  
```cd rsFQC```  
```cargo build --release```  
The resulting binary will be located in ```./rsFQC/target/release/rsFQC```  
You might need to flag the file as executable with ```chmod +x ./rsFQC```  

### Things to consider  

To speed things up, rsFQC will take shortcuts:  
- By testing only the first 50 nucleotides to get duplication levels, which might inflate duplication levels
- In long read mode, the quality chart is binned every 10 nucleotides.

### Usage
Only one argument needed - (no flags needed)

If you need to analyse multiple files at once in the same directory (multi mode):  
```rsFQC /path/to/fatqs/*```  
This will create a summary file in the current working directory named ```rsFQC.summary.txt```

If you need to analyse only one file:  
```rsFQC /path/to/file.fq.gz```  
This will create a in-terminal report, see below.


### Output example (single mode)
```rsFQC file.fq.gz```

Will produce an output that looks like this (piped to stdout so it's easy to capture in a file with ```rsFQC file.fq.gz > report.txt```:  
```
Detected gzip encoding
~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~       rsFQC       ~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sampling the first 100,000 records
of file
file.fastq.gz
-----------------------------------------------------
Long reads mode
-----------------------------------------------------
QUALITY
-----------------------------------------------------

y = Mean quality score at each position (horizontal line = Q20)
⡁                                                                  ⡇⢸⢸  90.0
⠄                                                                  ⡇⢸⢸ 
⠂                                                                  ⡇⢸⢸ 
⡁                                                                  ⡇⢸⢸ 
⠄                                                                  ⡇⢸⢸ 
⠂                                                                  ⡇⢸⢸ 
⡁                                                         ⢸     ⡇⡆ ⡇⢸⢸ 
⠄                                                     ⢀   ⢸⢰   ⢠⡇⡇ ⡇⢸⢸ 
⠂                                                     ⢸   ⢸⢸   ⢸⡇⡇ ⡇⢸⢸ 
⡁                                                ⡀⢰ ⢰ ⢸   ⣼⢸   ⢸⡇⡇⢀⣷⢸⢸⡀
⠄                                             ⣰ ⡀⣧⢸ ⢸⢀⣼ ⢀⡀⣿⣼⡀  ⢸⣷⣇⣾⣿⣸⣾⡇
⡂                                   ⡀⡀  ⢀ ⣀⡆⣠⣠⣿⣶⣇⣿⢸⣀⣼⣸⣿⣿⣾⣿⣿⣿⣿⣰⣀⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣭⣿⣿⣯⣭⣯⣿⣽⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇
⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇
⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠁ 3.0
0.0                                                                47621.0


y = Distribution of mean read quality
                 ⡠⠤⠒⡇                                                   18759.0
                ⡰⠁  ⠸⡀                                                 
               ⡰⠁    ⢣                                                 
              ⢠⠃     ⠈⡆                                                
             ⢠⠃       ⢸                                                
            ⢠⠃        ⠘⡄                                               
           ⢠⠊          ⢇                                               
          ⢠⠃           ⢸                                               
         ⢠⠃             ⡇                                              
        ⢠⠃              ⢱                                              
       ⢠⠃               ⠘⡄                                             
      ⢠⠃                 ⢱                                             
     ⢠⠃                   ⢇                                            
    ⡠⠃                    ⠘⡄                                           
  ⢠⠊                       ⠈⠑⢄                                         
⠉⠉⠁                           ⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠁ 1.0
11.0                                                               39.0

-----------------------------------------------------
Length
-----------------------------------------------------
Mean	Median	Min	Max
1693	6607	200	47621

y = Distribution of read length at each position
⡇                                                                       125.0
⡇                                                                      
⡇                                                                      
⣷                                                                      
⣿                                                                      
⣿                                                                      
⣿⡆                                                                     
⣿⡇                                                                     
⣿⡇                                                                     
⣿⣧⡀                                                                    
⣿⣿⡇                                                                    
⣿⣿⣇                                                                    
⣿⣿⣿⣦                                                                   
⣿⣿⣿⣿⣧⣤⡀                                                                
⣿⣿⣿⣿⣿⣿⣿⣶⣤⣤⣄⣠⣀⣀⣀⡀⢀⡀⡀                                                    
⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠁ 1.0
200.0                                                              47621.0

-----------------------------------------------------
Duplications
-----------------------------------------------------
Before dedup: 100,000 reads
After dedup: 96,729 reads
Reads are 3.27% duplicated 
96.73% left if deduplicated
y = Number of reads duplicated x times
⡇                                                                       1122.0
⡇                                                                      
⢇                                                                      
⢸                                                                      
⢸                                                                      
⠘⡄                                                                     
 ⡇                                                                     
 ⡇                                                                     
 ⢱                                                                     
 ⢸                                                                     
 ⢸                                                                     
  ⡇                                                                    
  ⢸                                                                    
   ⢇⡀                                                                  
    ⠈⠑⠢⢄⣀⣀⣀⣀⣀⣀⡀                                                        
              ⠈⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠁ 1.0
2.0                                                                40.0

```
