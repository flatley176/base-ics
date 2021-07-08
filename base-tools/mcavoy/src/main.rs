//use std::fs;
use clap::{Arg,App};
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

/*struct McFasta {
    records: Vec<McFastaRecord>
}*/

struct McFastaRecord {
    identifier: String,
    description: String,
    sequence: String
}

impl McFastaRecord {
    
    /*fn new_fasta_record(id: String, desc: String, seq: String) -> McFastaRecord {
        McFastaRecord {
            identifier: id,
            description: desc,
            sequence: seq
        }
    }*/
    
    fn new_fasta(file_path: &str) -> Vec<McFastaRecord> {
        let file_handle: File = File::open(file_path).unwrap();
        let buffered_file: BufReader<File> = BufReader::new(file_handle);
        let mut fasta_records = Vec::new();
        let mut current_rec;
        let mut current_rec_identifier: String = String::from("");
        let mut current_rec_description: String = String::from("");
        let mut current_rec_sequence: String = String::from("");
        for line in buffered_file.lines() {
            let entry = &mut line.unwrap();
            *entry = entry.trim_end().to_string();
            match entry.chars().next() {
                Some('>') => { // header line 
                                entry.remove(0);
                                match entry.contains(char::is_whitespace) {
                                true => {   // header line with white space
                                            println!("Start of new fasta sequence");
                                            //println!("Line is {}", entry);
                                            let header_list: Vec<&str> = entry.split_whitespace().collect();  
                                            println!("ID is {}, Description is {}", header_list[0], header_list[1]); 
                                            match fasta_records.len() { // empty list
                                                0 => {
                                                        println!("Reading sequences");
                                                     }
                                                _ => { // nonempty list
                                                        current_rec = McFastaRecord {
                                                                        identifier: current_rec_identifier,
                                                                        description: current_rec_description,
                                                                        sequence: current_rec_sequence                           
                                                        };
                                                        fasta_records.insert(fasta_records.len(), current_rec);
                                                        println!("Inserted {}", &current_rec.sequence);
                                                     }
                                            }
                                            current_rec_identifier = header_list[0].to_string();
                                            current_rec_description = header_list[1].to_string();
                                        }
                                  _  => { // header line, but no whitespace
                                            match fasta_records.len() {
                                                0 => { // empty list
                                                        println!("Reading sequences");
                                                     }                                            
                                                _ => { // records vector is not empty
                                                        current_rec = McFastaRecord { identifier: current_rec_identifier,
                                                                                      description: current_rec_description,
                                                                                      sequence: current_rec_sequence 
                                                        };
                                                        fasta_records.insert(fasta_records.len(), current_rec);
                                                        println!("Inserted {}", &current_rec.sequence);
                                                     }
                                            }
                                            current_rec_identifier = entry.to_string();
                                            current_rec_description = String::from("");
                                            println!("ID is {}", &entry);
                                        }
                                }
                                current_rec_sequence = String::from("");
                              }
                               
                       _  =>  { // sequence line
                                current_rec_sequence.push_str(entry);
                                //println!("Line is {}", entry);
                              }
                
            } 
        }
        fasta_records 
    }

    /*fn print_fasta(&self) {

    }

   fn get_length(&self) -> Vec<(String, u32)> {

    }
    
    fn get_gc_content(&self) -> Vec<(String, f32> {

    }

    fn get_base_composition(&self) -> Vec<(String, u32, u32, u32, u32)> {

    }
    
    fn get_sequence_complexity(&self, n_bases) -> Vec<(String, f32)> {
        
    }

    // moonshot: some kind of low impact BLAST that tells you what's in the fasta
    // minimally returns: query, subject, % identity and bitscore 
    fn whats_in_the_pot(&self) -> Vec<(String, String, f32, f32)> {

    }*/
} 



fn main() {
    let argument_set: clap::ArgMatches = App::new("AOC: Day 2")
                            .arg(Arg::with_name("input")
                                .short("f")
                                .long("fasta")
                                .takes_value(true)
                                .required(true))
                            .get_matches();
    let input_file: &str = argument_set.value_of("input").unwrap();
    let records = McFastaRecord::new_fasta(input_file);
}
