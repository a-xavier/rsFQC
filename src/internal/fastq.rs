use crate::internal::func::{*};
use rayon::prelude::*;
use itertools::Itertools;
use textplots::{Chart, Plot, Shape, ColorPlot};
use rgb::RGB8;

#[derive(Debug,Clone)]
pub struct FqObject  {
    // General
    pub filepath: String,
    pub number_of_records_used: usize,
    pub first_n_sequences: Vec<String>,
    pub first_n_quality: Vec<Vec<u32>>,
    pub gzipped: bool, // default is false before being processed
    pub processed: bool, // default is false before being processed
    pub isFastq: bool, // default is false before being processed
    pub isFile: bool,  // default is false before being processed
    pub isReadable: bool,  // default is false before being processed

    // Type of sequencing
    pub isLongReads: bool,


    // plotting
    pub plot_width: u32,
    pub plot_height: u32,

    // Quality
    pub qual_min: u32,
    pub qual_mean: u32,
    pub qual_median: u32,
    pub qual_max: u32,
    // Length
    pub len_min: u32,
    pub len_mean: u32,
    pub len_median: u32,
    pub len_max: u32,
    // Duplication
    pub duplication_levels: f32
}

impl FqObject{
    // Create new unprocessed 
    pub fn new(filepath: String, number_of_records_to_get: usize) -> Self {
        Self {
            filepath: filepath,
            isFile: false,
            number_of_records_used: number_of_records_to_get,
            first_n_sequences: Vec::new(),
            first_n_quality: Vec::new(),
            gzipped: false, 
            processed: false,
            isFastq: false,
            isReadable:false, 

            // Type of sequencing
            isLongReads: false,

            // plotting
            plot_width: 140,
            plot_height: 60,
            
            // Yet unknown
            // Quality
            qual_min: 0,
            qual_mean: 0,
            qual_median: 0,
            qual_max: 0,
            // Length
            len_min: 0,
            len_mean: 0,
            len_median: 0,
            len_max: 0,
            // Duplication
            duplication_levels: 0.0
        }
    }

    // Process the file
    pub fn pre_process(& mut self){

        // 0 is FILE - not directory
        self.isFile = check_file(&self.filepath);

        if self.isFile{
            // 1 - Get encoding and set flag - GZIP OR OTHER
            self.gzipped  = get_encoding(&self.filepath);

            // 2 - Can get read as UTF8
            self.isReadable = is_readable(&self.filepath, self.gzipped);

            if self.isFile & self.isReadable {

                //3  - Is it a proper FastQ file
                self.isFastq = is_fastq_file(&self);
            }
        }
    } // used to filter our fastq holder and remove bad files

    pub fn process_single(& mut self){
        // 3 populate the main fields
        self.first_n_sequences = get_first_n_sequence_lines_of_fastq_file(&self);
        self.first_n_quality = get_first_n_quality_lines_of_fastq_file(&self);
        
        self.length_quartiles();
        if self.len_max >= 1000 {self.isLongReads = true}; // set the longreads tag if long reads are detected
        header(self);
        sep();
        if self.isLongReads {println!("Long Reads Mode");}else{println!("Short Read Mode")}
        self.quality_quartiles();
        self.duplication_calculation();
        sep();
        println!("DUPLICATION");
        sep();
        println!("Duplication level: {}%", 100.-self.duplication_levels);
        self.duplication_chart();
        sep();
        println!("QUALITY");
        sep();
        println!("Mean Read Quality Distribution");
        println!("Min Q\tMed Q\tAvg Q\tMax Q");
        println!("{}\t{}\t{}\t{}", self.qual_min, self.qual_median, self.qual_mean, self.qual_max);
        self.quality_charts();
        sep();
        println!("LENGTH");
        sep();
        println!("Read Length Distribution");
        println!("Min L\tMed L\tAvg L\tMax L");
        println!("{}\t{}\t{}\t{}", self.len_min, self.len_median, self.len_mean, self.len_max);
        self.length_charts();
    }

    pub fn process_multi(& mut self){
        // 3 populate the main fields
        self.first_n_sequences = get_first_n_sequence_lines_of_fastq_file(&self);
        self.first_n_quality = get_first_n_quality_lines_of_fastq_file(&self);
        self.length_quartiles();
        if self.len_max >= 1000 {self.isLongReads = true}; // set the longreads tag if long reads are detected
        self.quality_quartiles();
        self.duplication_calculation();
    }

    pub fn length_quartiles(& mut self){

        let size_of_filtered_lines: Vec<usize> = self.first_n_sequences
        .par_iter()
        .map(|x: &String| x.len())
        .collect();

        let maximum = size_of_filtered_lines.par_iter().max().unwrap().to_owned();
        let minimum = size_of_filtered_lines.par_iter().min().unwrap().to_owned();

        let total : usize = size_of_filtered_lines.par_iter().sum();
        let mean = total / size_of_filtered_lines.len();
        let median = size_of_filtered_lines.get(size_of_filtered_lines.len()/2).unwrap().to_owned();

        self.len_min = minimum as u32;
        self.len_max = maximum as u32;
        self.len_median = median as u32;
        self.len_mean = mean as u32;
    }

    // QUality Distributions
    pub fn quality_quartiles(& mut self){

        let mean_quality_of_filtered_lines: Vec<u32> = self.first_n_quality
        .par_iter()
        .map(|x|  x.par_iter().sum::<u32>()/ (x.len() as u32))
        .collect();

        let maximum = mean_quality_of_filtered_lines.par_iter().max().unwrap().to_owned();
        let minimum = mean_quality_of_filtered_lines.par_iter().min().unwrap().to_owned();

        let total : u32 = mean_quality_of_filtered_lines.par_iter().sum();
        let mean = total / (mean_quality_of_filtered_lines.len() as u32);
        let median = mean_quality_of_filtered_lines.get(mean_quality_of_filtered_lines.len()/2).unwrap().to_owned();

        self.qual_min = minimum as u32;
        self.qual_max = maximum as u32;
        self.qual_median = median as u32;
        self.qual_mean = mean as u32;
    }

    pub fn duplication_calculation(& mut self) {
                // Limit to 50 bp? like fastqc
                let filtered_lines_to_get: Vec<String> = self.first_n_sequences
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
                
                // println!("Before dedup: {} reads", pretty_print_int(&full_length));
                // println!("After dedup: {} reads", pretty_print_int(&dedup_length));
            
                // println!("Reads are {:02.2}% duplicated ", 100.0 - dedup_percent);
                // println!("{:02.2}% left if deduplicated", dedup_percent);

                self.duplication_levels = dedup_percent;
    }

    // DIsplay a chart for duplication
    pub fn duplication_chart(& mut self) -> () {
                 // Limit to 50 bp? like fastqc
                 let filtered_lines_to_get: Vec<String> = self.first_n_sequences
                 .par_iter()
                 .map(|x: &String| 
                     if x.len() < 50{
                         x.to_owned()}else{
                         x[0..50].to_owned()}
                 )
                 .collect();

                // itertools test
                // Gives the number of times a sequence appear in our sample
                // {(2, "ATGCGTCG"), 
                //  (4, "TCGCCGCGGA"),
                //  (1, "TCGCGCGAGGNTTA")}
                let mut vector_for_map = filtered_lines_to_get.clone();
                vector_for_map.sort();
                let count_map: Vec<(usize, String)> = vector_for_map
                .into_iter()
                .dedup_with_count()
                .collect();
            
                // Todo -> Bin all reads duplicated more than 10 times
                // ONLY KEEP TUPLES WHEN ELEMENT 0 is more than 1 - Duplicated sequences
                let dup_map: Vec<(usize, String)> = count_map
                .into_par_iter()
                .filter(|x| x.0 > 1)
                .collect();
                

                // Only keep the first element of each tuple 
                // Result should be {1, 2, 3, 2, 3, 4, 4, 1, ,1 ,1, 4, 2}
                // Replace anything over 10 by 10
                let occurence_vec: Vec<usize> = dup_map
                .par_iter()
                .map(|(number, _)| if number > &10 {10 as usize} else {number.to_owned()} )
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
                    Chart::new(self.plot_width, self.plot_height, 2.0, max_occurrence as f32) // Start duplication chart at x = 2 
                    .lineplot(&Shape::Lines(&points_for_distribution_chart))
                    .display();
                }
        }

        // Display Quality Charts
    pub fn quality_charts(& mut self) -> () {
        let mut step = 1;
            if self.isLongReads {step = 10} // if long read, sample every 10 base
            
            // Number of reads with that read length
            let population_numbers: Vec<u32> = (1..self.len_max +1).step_by(step)
            .into_iter()
            .map(|x| sub_vector_above_this_length(&self.first_n_quality, x as usize))
            .collect();
        
            // Sum of qualities for each positions
            let sum_numbers: Vec<u32> = (0..self.len_max).step_by(step)
            .into_iter()
            .map(|x| sub_vector_sum_at_index(&self.first_n_quality, x as usize))
            .collect();

            // Get the points for chart
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
            Chart::new(self.plot_width, self.plot_height, 0.0, self.len_max as f32)
            .lineplot(&Shape::Bars(&quality_points_for_chart))
            .linecolorplot(&Shape::Continuous(Box::new(|_x| qual_threshold)), red)
            .display();
        
            // 2nd chart with mean quality per read
            let mut mean_quality_vector: Vec<u32>  = self.first_n_quality
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
            Chart::new(self.plot_width, self.plot_height, min_quality as f32, max_quality as f32)
            .lineplot(&Shape::Lines(&point_for_mean_qual))
            .display();
        }

        pub fn length_charts(& mut self) -> () {
            let size_of_filtered_lines: Vec<usize> = self.first_n_sequences
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
            Chart::new(self.plot_width, self.plot_height, self.len_min as f32, self.len_max as f32)
                .lineplot(&Shape::Bars(&distribution_points))
                .display();

        }

}
