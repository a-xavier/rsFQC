#![allow(non_snake_case)]
use std::{env, usize};
// Inside crates
mod internal;
use crate::internal::fastq::FqObject;
use crate::internal::func::{*};

fn main() {

    // DEBUG
    env::set_var("RUST_BACKTRACE", "1");

    // Input
    let mut args: Vec<String> = env::args().collect();

    // Capture the executable location
    let _executable_location = args.get(0).unwrap().to_owned();

    // Test if there is any argument used
    args.remove(0); // remove first argument since it's the location of the executable 
    if args.len() == 0 {
        println!("No input");
        std::process::exit(1);
    };

    // The number of records to test
    let number_of_records_to_test: usize = 100000;

    // CHART WIDTH
    let _chart_width = 140 as u32;
    let _chart_height = 60 as u32;

    //single of multi mode
    let mut multi_mode: bool = false;

    // Get the list of Files to process
    let mut all_fq_to_process: Vec<FqObject> = Vec::new();
    for path in args {
        let mut new_fq = FqObject::new(path, number_of_records_to_test);
        new_fq.pre_process();
        if new_fq.isFastq & new_fq.isFile {all_fq_to_process.push(new_fq)}
    }

    // If no valid FastQ detected
    if all_fq_to_process.len() == 0 {
        println!("No valid input detected.");
        std::process::exit(1);
    };

    // Check the number of Fastq file to process
    if all_fq_to_process.len() > 1 {multi_mode = true;println!("Found {} valid FastQ files.", all_fq_to_process.len())};

    // Process all Fastq Files
    let mut new_holder: Vec<FqObject> = Vec::new();
    if multi_mode {
        println!("Processing.");
        for mut fq in  all_fq_to_process {
            fq.process_multi();
            new_holder.push(fq)}
            println!("Processed {} FastQ files into {}/rsFQC.summary.txt", new_holder.len(), env::current_dir().unwrap().to_str().unwrap());
        write_reports(new_holder);
    } else{
        all_fq_to_process.get(0).unwrap().to_owned().process_single();
    }


    // FASTQ FORMAT
    // @SEQ_ID
    // GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    // +
    // !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

    // // EXIT CODES

    std::process::exit(0);

}