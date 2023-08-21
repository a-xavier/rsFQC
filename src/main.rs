#![allow(non_snake_case)]
use std::fs::File;
use std::io::{BufReader, BufRead, Read};
use std::{env, usize};
// Outside crates
use flate2::read::GzDecoder;
use rayon::prelude::*;
use textplots::{Chart, Plot, Shape, ColorPlot};
use rgb::RGB8;
use itertools::Itertools;

fn main() {
    // Input
    let args: Vec<String> = env::args().collect();
    // 0 is the name of the binary
    let filepath = args.get(1).unwrap_or_else(||{ 
        println!("No input");
        std::process::exit(1);
    });

    // The number of records to test
    let number_of_records_to_test: usize = 100000;

    // Iterate over lines and process them
    // get first nth lines 
    
    // Get indices to get.
    // second lines of group of 4
    // if zero based: 1 - 5 - 9 - 13 - etc
    // FASTQ FORMAT
    // @SEQ_ID
    // GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    // +
    // !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

    // PRINT STUFF
    let lines_to_test: Vec<String> = buffer_to_fq_lines(filepath, number_of_records_to_test);
    
    // Check if respect fastq format
    let first = lines_to_test[0].chars().nth(0).unwrap();
    let second = lines_to_test[2].chars().nth(0).unwrap();

    // Exit if not fastq file
    if first != '@' && second != '+' {
        println!("This is not a valid FastQ file. (We test for '@' on first line and '+' on third line)");
        std::process::exit(1);
    }

    let (mean, median, maximum, minimum) = mean_and_median(&lines_to_test, number_of_records_to_test);

    // Check if long reads or short reads
    let mut length_mode: &str= "Short reads";
    if maximum > 1000 {
        length_mode = "Long reads"
    }

    header(&number_of_records_to_test, filepath);
    sep();
    println!("{} mode", length_mode);
    
    //  QUALITY
    sep();
    println!("QUALITY");
    sep();
    quality_chart(&lines_to_test, number_of_records_to_test, &maximum, &minimum, length_mode);
    
    //  LENGTH
    sep();
    println!("Length");
    sep();
    println!("Mean\tMedian\tMin\tMax");
    println!("{}\t{}\t{}\t{}", mean, median, minimum, maximum);
    length_distribution_chart(&lines_to_test, number_of_records_to_test, &maximum, &minimum);
    // Duplication
    sep();
    println!("Duplications");
    sep();
    calculate_duplication(&lines_to_test, number_of_records_to_test, &maximum, &minimum);

    // EXIT CODES

    std::process::exit(0);

}

/// Function to detect the encoding of fq file -> returns either "gzip" or "text"
fn get_encoding(filepath: &String) -> &str{
    //Use first 2 bytes (byte1 == 0x1f) && (byte2 == 0x8b) 
    let mut f = File::open(filepath).unwrap();
    let mut buf = [0u8; 2];
    f.read_exact(&mut buf).unwrap();

    let byte_1 = format!("0x{:02x}", buf[0]);
    let byte_2 = format!("0x{:02x}", buf[1]);

    if (byte_1 == "0x1f") && (byte_2 == "0x8b"){
        println!("Detected gzip encoding");
        return "gzip"
    } else{
        println!("Detected plain text encoding");
        return "text"
    }
}

/// Given a filename and a number of records to get, returns a Vec<String> of fastq lines
fn buffer_to_fq_lines(filepath: &String, number_of_records_to_get: usize)-> Vec<String>{
    // Open file in any case
    let file = File::open(filepath).unwrap();
    // Initialise the holder of lines to test
    let mut lines_to_test:Vec<String> = Vec::new(); 
    // Get proper encoding
    // TODO: DO NOT DO TWO CONNECTIONS /OPEN 1 to test encoding, one to read file
    let encoding = get_encoding(filepath);
    if encoding == "gzip"{
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

/// get mean and median
fn mean_and_median(lines_to_test: &Vec<String>, number_of_records_to_get: usize) -> (usize, usize, usize, usize){
    let index_to_get: Vec<usize> = (1..(number_of_records_to_get*4)).step_by(4).collect();

    let filtered_lines_to_get: Vec<String> = index_to_get
    .par_iter()
    .filter_map(|&index| lines_to_test.get(index).cloned())
    .collect();

    let size_of_filtered_lines: Vec<usize> = filtered_lines_to_get
    .par_iter()
    .map(|x: &String| x.len())
    .collect();

    let maximum = size_of_filtered_lines.par_iter().max().unwrap().to_owned();
    let minimum = size_of_filtered_lines.par_iter().min().unwrap().to_owned();

    let total : usize = size_of_filtered_lines.par_iter().sum();
    let mean = total / size_of_filtered_lines.len();
    let median = size_of_filtered_lines.get(size_of_filtered_lines.len()/2).unwrap().to_owned();
    // let formatted_mean = format!("{}/{}\n", mean, median);

    return (mean, median, maximum, minimum);
}

///quality chart
fn quality_chart(lines_to_test: &Vec<String>, number_of_records_to_get: usize, maximum: &usize, _minimum: &usize ,lenght_mode: &str){
    let index_to_get: Vec<usize> = (3..(number_of_records_to_get*4)).step_by(4).collect();
    
    let filtered_lines_to_get: Vec<String> = index_to_get.par_iter().filter_map(|&index| lines_to_test.get(index).cloned()).collect();

    let mut step = 1;

    if lenght_mode == "Long reads" {step = 10} // if long read, sample every 10 base
    
    // For each line of zplkodjdjqwe quality , get the vector of qualities [ [20, 32, 43], [32, 43, 13], [12, 32, 43] ]
    let quality_int_lines: Vec<Vec<u32>>  = filtered_lines_to_get
    .par_iter()
    .map(|x: &String| quality_vector_from_line(x.to_owned()))
    .collect();
    
    // Number of reads with that read length
    let population_numbers: Vec<u32> = (1..*maximum+1).step_by(step)
    .into_iter()
    .map(|x| sub_vector_above_this_length(&quality_int_lines, x))
    .collect();

    // Sum of qualities for each positions
    let sum_numbers: Vec<u32> = (0..*maximum).step_by(step)
    .into_iter()
    .map(|x| sub_vector_sum_at_index(&quality_int_lines, x))
    .collect();

    let quality_points_for_chart: Vec<(f32, f32)> = (0..sum_numbers.len())
    .into_iter()
    .map(|x| (x as f32 * step as f32 + 1.0, (sum_numbers.get(x).unwrap().to_owned() as f32 / population_numbers.get(x).unwrap().to_owned() as f32))).
    collect();

    // Red
    let red = RGB8 {r:100, g:255, b:0};
    let qual_threshold = 20 as f32;

    println!("\ny = Mean quality score at each position (horizontal line = Q20)");

    // &Shape::Lines(&[(0.0, qual_threshold), (maximum.to_owned() as f32, qual_threshold)])
    // 
    Chart::new(180, 60, 0.0, maximum.to_owned() as f32)
    .lineplot(&Shape::Bars(&quality_points_for_chart))
    .linecolorplot(&Shape::Continuous(Box::new(|_x| qual_threshold)), red)
    .display();

    // 2nd chart with mean quality per read
    let mut mean_quality_vector: Vec<u32>  = quality_int_lines
    .par_iter()
    .map(|x| x.par_iter().sum::<u32>() / x.len() as u32)
    .collect();
    mean_quality_vector.sort();

    let vec_for_extremes = mean_quality_vector.clone(); // cant get around that borrowed thing

    let table_dup: Vec<(usize, u32)> = mean_quality_vector.into_iter().dedup_with_count().collect();
    let point_for_mean_qual: Vec<(f32, f32)> = table_dup
    .into_par_iter()
    .map(|(x, y)| (y as f32, x as f32))
    .collect();

    let max_quality: u32 = vec_for_extremes.iter().cloned().max().unwrap();

    let min_quality: u32 = vec_for_extremes.iter().cloned().min().unwrap();

    // dbg!(&point_for_mean_qual);

    println!("\ny = Distribution of mean read quality");
    // &Shape::Lines(&[(0.0, qual_threshold), (maximum.to_owned() as f32, qual_threshold)])
    Chart::new(180, 60, min_quality as f32, max_quality as f32)
    .lineplot(&Shape::Lines(&point_for_mean_qual))
    .display();


}

/// Character to quality
fn char_to_qual(c: char)-> u32{
    let ascii_code = c as u32;
    return ascii_code-33
}

/// Quality line to vector of quality
fn quality_vector_from_line(line: String)-> Vec<u32>{
    let return_value: Vec<u32> = line.chars().into_iter().par_bridge().map(|x| char_to_qual(x)).collect();
    return return_value
}

/// Get the number of sub vectors above a certain length
fn sub_vector_above_this_length(initial_vec: &Vec<Vec<u32>>, length_to_test: usize) -> u32{
    let caca: usize = initial_vec.par_iter().map(|x| x.len() >= length_to_test).collect::<Vec<bool>>()
    .par_iter().filter(|x| **x).count();
    return caca as u32;
}

/// Get sum of subvectors at index x
fn sub_vector_sum_at_index(initial_vec: &Vec<Vec<u32>>, index: usize) -> u32{
    let caca: Vec<u32> = initial_vec
    .par_iter()
    .map(|x| x.get(index).unwrap_or(&0).to_owned())
    .collect();
    let sumsum: u32 = caca.par_iter().sum();
    return sumsum
}

/// Length distribution
fn length_distribution_chart(initial_vec: &Vec<String>, number_of_records_to_get: usize, maximum: &usize, minimum: &usize){
    // GET ATGC LINE // Line 2 if 1 based /
    let index_to_get: Vec<usize> = (1..(number_of_records_to_get*4)).step_by(4).collect();

    let filtered_lines_to_get: Vec<String> = index_to_get
    .par_iter()
    .filter_map(|&index| initial_vec.get(index).cloned())
    .collect();

    let size_of_filtered_lines: Vec<usize> = filtered_lines_to_get
    .par_iter()
    .map(|x: &String| x.len())
    .collect();   

    let mut unique_lengths = size_of_filtered_lines.clone();
    unique_lengths.sort();
    unique_lengths.dedup();

    let distribution_points: Vec<(f32, f32)> = unique_lengths
    .par_iter()
    .map(|x| (*x as f32, count_number_of_occurence(&size_of_filtered_lines, x.to_owned())))
    .collect();

    println!("\ny = Distribution of read length at each position");
    Chart::new(180, 60, minimum.to_owned() as f32, maximum.to_owned() as f32)
        .lineplot(&Shape::Bars(&distribution_points))
        .display();


}

/// Count number of occurence of something in a vector
fn count_number_of_occurence(initial_vec: &Vec<usize>, item: usize)-> f32{
    let tmp: usize = initial_vec
    .par_iter()
    .map(|x: &usize| x.to_owned() == item)
    .collect::<Vec<bool>>().par_iter().filter(|x: &&bool| **x).count();

    return tmp as f32;
}

/// Show header
fn header(number_of_records_to_get: &usize, filename: &String){
    println!("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    println!("~~~~       rsFQC       ~~~~");
    println!("~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    println!("Sampling the first {} records", pretty_print_int(number_of_records_to_get));
    println!("of file");
    println!("{}",filename);
}

/// print separator
fn sep(){
    println!("-----------------------------------------------------")
}

/// Duplications (Chart?)

fn calculate_duplication(initial_vec: &Vec<String>, number_of_records_to_get: usize, _maximum: &usize, _minimum: &usize){
    let index_to_get: Vec<usize> = (1..(number_of_records_to_get*4)).step_by(4).collect();


    let mut filtered_lines_to_get: Vec<String> = index_to_get
    .par_iter()
    .filter_map(|&index| initial_vec.get(index).cloned())
    .collect();

    // Limit to 50 bp? like fastqc
    filtered_lines_to_get = filtered_lines_to_get
    .par_iter()
    .map(|x: &String| 
        if x.len() < 50{
            x.to_owned()}else{
            x[0..50].to_owned()}
    )
    .collect();

    let mut unique_lines = filtered_lines_to_get.clone();
    unique_lines.sort();
    unique_lines.dedup();

    let full_length = filtered_lines_to_get.len();
    let dedup_length = unique_lines.len();

    // Duplication levels
    let dedup_percent: f32 = (dedup_length as f32 /full_length as f32) * 100.0;

    println!("Before dedup: {} reads", pretty_print_int(&full_length));
    println!("After dedup: {} reads", pretty_print_int(&dedup_length));

    println!("Reads are {:02.2}% duplicated ", 100.0 - dedup_percent);
    println!("{:02.2}% left if deduplicated", dedup_percent);


    // itertools test
    let mut vector_for_map = filtered_lines_to_get.clone();
    vector_for_map.sort();
    let count_map: Vec<(usize, String)> = vector_for_map
    .into_iter()
    .dedup_with_count()
    .collect();

    // Todo -> Bin all reads duplicated more than 10 times
    let dup_map: Vec<(usize, String)> = count_map
    .into_par_iter()
    .filter(|x| x.0 > 1)
    .collect();

    let occurence_vec: Vec<usize> = dup_map
    .par_iter()
    .map(|(number, _)| number.to_owned() )
    .collect();

    let mut unique_occurrences: Vec<usize> = occurence_vec.clone();
    unique_occurrences.sort();
    unique_occurrences.dedup();

    
    if unique_occurrences.is_empty() {
        println!("No duplication detected!")
    } else {

    // PANIC IF NO DUPLICATED AT ALL
        let max_occurrence = unique_occurrences.iter().max().unwrap().to_owned();
        let points_for_distribution_chart: Vec<(f32, f32)> = unique_occurrences
        .par_iter()
        .map(|x| (*x as f32, count_number_of_occurence(&occurence_vec, *x) as f32))
        .collect();


        println!("y = Number of reads duplicated x times");
        Chart::new(180, 60, 2.0, max_occurrence as f32) // Start duplication chart at x = 2 
        .lineplot(&Shape::Lines(&points_for_distribution_chart))
        .display();
    }
    



    // Brute force because i don't know how to do it efficiently
    // let mut unique_lines = filtered_lines_to_get.clone();
    // unique_lines.sort();
    // unique_lines.dedup();
    // let dup_vector_1 = count_number_of_occurence_string_in_other_vector(&unique_lines, &filtered_lines_to_get);
    // dbg!(dup_vector_1);
}

/// Count number of time each element of a vector appear in another vector
/// Brute force because i don't know how to do it efficiently
fn _count_number_of_occurence_string_in_other_vector(initial_vec: &Vec<String>, other_vector: &Vec<String>)-> Vec<usize>{
    let tmp: Vec<usize> = initial_vec
    .par_iter()
    .map(|x: &String| 
        other_vector.par_iter().filter(|s| s.to_owned() == x).count()
    )
    .collect::<Vec<usize>>();

    return tmp;
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