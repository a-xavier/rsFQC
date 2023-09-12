use std::{fs::File, io::{Read, BufReader, BufRead, Write}, path::PathBuf};
use flate2::read::GzDecoder;
use super::fastq::FqObject;
use rayon::prelude::*;

/// Function to detect the encoding of fq file -> returns either "gzip" or "text"
pub fn get_encoding(filepath: &String) -> bool{
    //Use first 2 bytes (byte1 == 0x1f) && (byte2 == 0x8b) 
    let mut f = File::open(filepath).unwrap();
    let mut buf = [0u8; 2];
    f.read_exact(&mut buf).unwrap();

    let byte_1 = format!("0x{:02x}", buf[0]);
    let byte_2 = format!("0x{:02x}", buf[1]);

    if (byte_1 == "0x1f") && (byte_2 == "0x8b"){
        // println!("Detected gzip encoding");
        return true
    } else{
        // println!("Detected plain text encoding");
        return false
    }
}

/// Function to test quickly if the file is actually a fastQ file
pub fn is_fastq_file(fqobject: &FqObject) -> bool {
    // Get first 3 lines
    let lines_to_test = buffer_to_fq_lines(&fqobject.filepath, 3, fqobject.gzipped);
    // Check if respect fastq format
    let first = lines_to_test[0].chars().nth(0).unwrap_or_else(|| '.');
    let second = lines_to_test[2].chars().nth(0).unwrap_or_else(|| '.');
    // Exit if not fastq file
    if first == '@' && second == '+' {
        return true;
    } else {
        return false;
    }
}

/// Given a filename , a number of records to get and a gzip flag
/// returns a Vec<String> of fastq lines
pub fn buffer_to_fq_lines(filepath: &String, number_of_records_to_get: usize, gzipped: bool)-> Vec<String>{
    // Open file in any case
    let file = File::open(filepath).unwrap();
    // Initialise the holder of lines to test
    let mut lines_to_test:Vec<String> = Vec::new(); 
    // Get proper encoding
    // TODO: DO NOT DO TWO CONNECTIONS /OPEN 1 to test encoding, one to read file
    if gzipped{
        let decoder: GzDecoder<File> = GzDecoder::new(file);
        // Create a buffered reader to read lines efficiently
        let reader: BufReader<GzDecoder<File>> = BufReader::new(decoder);
        lines_to_test.extend(reader.lines().take( number_of_records_to_get * 4 ).map(|x: Result<String, std::io::Error>| x.unwrap()).collect::<Vec<String>>());
    } else{
        let reader: BufReader<File> = BufReader::new(file);
        lines_to_test.extend(reader.lines().take( number_of_records_to_get * 4 ).map(|x: Result<String, std::io::Error>| x.unwrap()).collect::<Vec<String>>());
    } 
    return lines_to_test
}

/// Get the first n quality records of a fastq file
pub fn get_first_n_quality_lines_of_fastq_file(fqobject: &FqObject) -> Vec<Vec<u32>> {
    // Get first N records
    let lines_to_test: Vec<String> = buffer_to_fq_lines(&fqobject.filepath, fqobject.number_of_records_used, fqobject.gzipped);
    // Get the index of lines we need // Modulo 3 for quality lines
    let index_to_get: Vec<usize> = (3..(fqobject.number_of_records_used*4)).step_by(4).collect();  
    // Filter first N records based on the index of previous step
    let filtered_lines_to_get: Vec<String> = index_to_get.par_iter().filter_map(|&index| lines_to_test.get(index).cloned()).collect();

    let filtered_lines_to_quality: Vec<Vec<u32>> = filtered_lines_to_get.par_iter().map(|x| quality_vector_from_line(x.to_string())).collect();

    return filtered_lines_to_quality;
}

/// Get the first n sequence records of a fastq file
pub fn get_first_n_sequence_lines_of_fastq_file(fqobject: &FqObject) -> Vec<String>{
    // Get first N records
    let lines_to_test: Vec<String> = buffer_to_fq_lines(&fqobject.filepath, fqobject.number_of_records_used, fqobject.gzipped);
    // Get the index of lines we need // Modulo 1 for sequence lines
    let index_to_get: Vec<usize> = (1..(fqobject.number_of_records_used*4)).step_by(4).collect();  
    // Filter first N records based on the index of previous step
    let filtered_lines_to_get: Vec<String> = index_to_get.par_iter().filter_map(|&index| lines_to_test.get(index).cloned()).collect();

    return filtered_lines_to_get;
}

/// Character to quality
pub fn char_to_qual(c: char)-> u32{
    let ascii_code = c as u32;
    return ascii_code-33
}

/// Quality line to vector of quality
pub fn quality_vector_from_line(line: String)-> Vec<u32>{
    let return_value: Vec<u32> = line.chars().into_iter().par_bridge().map(|x| char_to_qual(x)).collect();
    return return_value
}

/// CHeck if it's a file and if I can read it
/// Is a bit clunky so far
pub fn check_file(filepath: &String) -> bool {
    let mut is_file = true;

    // Weed out directory
    if  PathBuf::from(filepath).is_dir() {is_file = false};

    return is_file

}

/// Check if It's readable into utf8
///  TODO: Find a better way
pub fn is_readable(filepath: &String, gzipped: bool) -> bool {
    let mut _is_it_readable: bool = true;
    if gzipped {
        let decoder: GzDecoder<File> = GzDecoder::new(File::open(filepath).unwrap());
        // Create a buffered reader to read lines efficiently
        let reader: BufReader<GzDecoder<File>> = BufReader::new(decoder);
        // Get 10th line and see if it's readable
        let result: Result<String, std::io::Error> = reader.lines().take(10).collect();
        _is_it_readable = match result {
            Ok(_file) => {true},
            Err(_e) => {false},
        }
    }else{
        // Weed out stuff that can't be read
        let reader: BufReader<File> = BufReader::new(File::open(filepath).unwrap());
        // Get 10th line and see if it's readable
        let result: Result<String, std::io::Error> = reader.lines().take(10).collect();
        _is_it_readable = match result {
            Ok(_file) => {true},
            Err(_e) => {false},
        }
    }   
    return _is_it_readable 
}

/// Stolen from Michael Hall https://stackoverflow.com/questions/26998485/is-it-possible-to-print-a-number-formatted-with-thousand-separator-in-rust
pub fn pretty_print_int(i: &usize) -> String {
    let mut s = String::new();
    let i_str = i.to_string();
    let a = i_str.chars().rev().enumerate();
    for (idx, val) in a {
        if idx != 0 && idx % 3 == 0 {
            s.insert(0, ',');
        }
        s.insert(0, val);
    }
    return s;
}

/// Count number of occurence of something in a vector
pub fn count_number_of_occurence(initial_vec: &Vec<usize>, item: usize)-> f32{
    let tmp: usize = initial_vec
    .par_iter()
    .map(|x: &usize| x.to_owned() == item)
    .collect::<Vec<bool>>().par_iter().filter(|x: &&bool| **x).count();
    return tmp as f32;
}

/// Get the number of sub vectors above a certain length
pub fn sub_vector_above_this_length(initial_vec: &Vec<Vec<u32>>, length_to_test: usize) -> u32{
    let caca: usize = initial_vec.par_iter().map(|x| x.len() >= length_to_test).collect::<Vec<bool>>()
    .par_iter().filter(|x| **x).count();
    return caca as u32;
}

/// Get sum of subvectors at index x
pub fn sub_vector_sum_at_index(initial_vec: &Vec<Vec<u32>>, index: usize) -> u32{
    let caca: Vec<u32> = initial_vec
    .par_iter()
    .map(|x| x.get(index).unwrap_or(&0).to_owned())
    .collect();
    let sumsum: u32 = caca.par_iter().sum();
    return sumsum
}

/// Write reports when in multi mode
pub fn write_reports(input: Vec<FqObject>) -> () {
    let mut file = File::create("rsFQC.summary.txt").unwrap();
    // file.write_all(b"Hello, world!").unwrap();
    file.write_all(b"File\tMinimum Length\tMedian Length\tAverage Length\tMaximum Lemgth\tMinimum Quality\tMedian Quality\tAverage Quality\tMaximum Quality\tDuplication Level\n").unwrap();
    for fq in input{
        let formated_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            fq.filepath,fq.len_min, fq.len_median, fq.len_mean, fq.len_max,
            fq.qual_min, fq.qual_median, fq.qual_mean, fq.qual_max, 
            fq.duplication_levels );   
        file.write_all(formated_line.as_bytes()).unwrap();
    }
}

/// Show header
pub fn header(fq:& mut  FqObject){
    println!("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    println!("~~~~       rsFQC       ~~~~");
    println!("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    println!("Sampling the first {} records", pretty_print_int(&fq.number_of_records_used));
    println!("of file");
    println!("{}",&fq.filepath);
}

/// print separator
pub fn sep(){
    println!("-----------------------------------------------------")
}