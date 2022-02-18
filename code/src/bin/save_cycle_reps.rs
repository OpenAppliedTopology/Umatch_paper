use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::clique::{Simplex, CliqueComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind, Indexing};

extern crate csv;

#[macro_use]
extern crate npy_derive;
extern crate npy;

use ndarray::{Array2, Array3};
use ndarray_npy::{read_npy, write_npy};
use tuple_conv::RepeatedTuple;

use std::io::{Read, Write, BufReader, BufRead};

use npy::NpyData;
use csv::ReaderBuilder;
use std::error::Error;
use std::fs::File;
use std::env;

use math::round;
use ordered_float::OrderedFloat;
use num::rational::Ratio;

use std::time::Instant;

fn main() -> Result<(), Box<dyn Error>> {
///////////////////////////////////////////////////////////////////////////////////////////////////
// Accessing commond line arguments
// Samle args:
// target\release\save_cycle_reps C:\Users\haibi\Documents\data_only_gdrive_upload\for_paper\clique\torus\dismat.npy 1 2 0.0 false
	let args: Vec<String> = env::args().collect();
	let mut data_file_name = args[1].clone();  // data file (.npy dissimilarity file) location
	let dim: usize = args[2].trim().parse().expect("Please type a number!"); // the dimension at which we want to compute the barcode
	let field: usize = args[3].trim().parse().expect("Please type a number!"); // the field parameter
	let mut maxdis: OrderedFloat<f64> = args[4].trim().parse().expect("Please type a number!"); // cutoff distance. 0.0 indicates radius
	let print_or_not: bool = args[5].trim().parse().expect("Please type a booliean!"); // print barcodes or not
	//let mut dir = args[6].clone(); // where to save the cycle representatives

	println!("Reading dissimilarity matrix data from: {:?}", data_file_name);
	let arr: Array2<f64> = read_npy(&data_file_name).unwrap();
/////////////////// Read dissimilarity matrix data into a Vec<Vec<FilVal>>

	let mut dis_mat = Vec::new();
	let mut min_max = OrderedFloat(0.0);
	for row in arr.outer_iter() {
		let mut vector = Vec::new();
		let mut max = OrderedFloat(0.0);
		for entry in row.iter() {
			let rounded_entry = round::floor(*entry, 15);
			vector.push(OrderedFloat(rounded_entry));
			if OrderedFloat(rounded_entry) > max { max = OrderedFloat(rounded_entry); }
		}
		dis_mat.push(vector);
		if max < min_max || min_max == OrderedFloat(0.0) { min_max = max; }
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

    let ringmetadata = RingMetadata{
        ringspec: RingSpec::Modulus(field),
        identity_additive: 0,
        identity_multiplicative: 1,
    };

	if maxdis == OrderedFloat(0.0) { maxdis = min_max; }
	let chx = CliqueComplex {
		dissimilarity_matrix: dis_mat,
		dissimilarity_value_max: maxdis,
		safe_homology_degrees_to_build_boundaries: (1..(dim+1)).collect(),
		major_dimension: MajorDimension::Row,
		ringmetadata,
		simplex_count: Vec::new()
	};

///////////////////////////////////////////////////////////////////////////////////////////////////

	println!("Caculating barcodes ...");
    let now = Instant::now();
    let blocks = factor_chain_complex(&chx, dim+1);
    println!("Runs {} seconds", now.elapsed().as_secs());
    println!("Done!");

    let mut barcode = blocks.barcode(dim);
    barcode.sort();
    println!("Barcode contains {} bars.", barcode.len());
    if print_or_not == true {
        for (start, end) in barcode.iter() {
            println!("{},{}", start, end);
        }
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

	use std::fs;

	let mut dir = String::from("reps_dim_");
	dir.push_str(&dim.to_string());
	fs::create_dir(dir.clone())?;
	println!("Saving cycle representatives in: {}", dir);

	let indexing = &blocks.dim_indexing[dim];
	let mut indexing_upper = &blocks.dim_indexing[dim+1];

	let mut counter = 0;
	for key in chx.keys_unordered_itr(dim) {
		if indexing.minkey_2_index.contains_key(&key) { continue; }
		else if indexing_upper.majkey_2_index.contains_key(&key) {
			let ind = indexing_upper.majkey_2_index[&key];
			let matched_key = &indexing_upper.index_2_minkey[ind];
			let diam1 = chx.key_2_filtration(&key);
			let diam2 = chx.key_2_filtration(matched_key);
			if diam1 == diam2 { continue; }

			let mut save_as = dir.clone();
			save_as.push_str("/rep_no_");
			save_as.push_str(&counter.to_string());
			save_as.push_str(".txt");
			let mut file = File::create(save_as).unwrap();

			let OrderedFloat(val) = key.filvalue;
			let mut tmp = serde_json::to_string(&key.vertices).unwrap();
			tmp += " ";
			tmp += &serde_json::to_string(&val).unwrap();
			tmp = String::from("Birth simplex/time: ") + &tmp + "\n";
			file.write(tmp.to_string().as_bytes()).unwrap();

			let OrderedFloat(val) = matched_key.filvalue;
			let mut tmp = serde_json::to_string(&matched_key.vertices).unwrap();
			tmp += " ";
			tmp += &serde_json::to_string(&val).unwrap();
			tmp = String::from("Death simplex/time: ") + &tmp + "\n";
			file.write(tmp.to_string().as_bytes()).unwrap();

			let mut tmp = String::from("Cycle representative combination:\nSimplicies Coefficients") + "\n";
			file.write(tmp.to_string().as_bytes()).unwrap();

			let hash = blocks.get_matched_basis_vector(dim, &key);
			for (simp, coef) in hash.iter(){
				let mut tmp = serde_json::to_string(&simp.vertices).unwrap();
				tmp += " ";
				tmp += &serde_json::to_string(&coef).unwrap();
				tmp += "\n";
				file.write(tmp.to_string().as_bytes()).unwrap();
			}
			counter += 1;
		} else {
			let mut save_as = dir.clone();
			save_as.push_str("/rep_no_");
			save_as.push_str(&counter.to_string());
			save_as.push_str(".txt");
			let mut file = File::create(save_as).unwrap();

			let OrderedFloat(val) = key.filvalue;
			let mut tmp = serde_json::to_string(&key.vertices).unwrap();
			tmp += " ";
			tmp += &serde_json::to_string(&val).unwrap();
			tmp = String::from("Birth simplex/time: ") + &tmp + "\n";
			file.write(tmp.to_string().as_bytes()).unwrap();

			let mut tmp = String::from("Death simplex/time: None\n");
			file.write(tmp.to_string().as_bytes()).unwrap();

			let mut tmp = String::from("Cycle representative:") + "\n";
			file.write(tmp.to_string().as_bytes()).unwrap();

			let hash = blocks.get_matched_basis_vector(dim, &key);
			for (simp, coef) in hash.iter(){
				let mut tmp = serde_json::to_string(&simp.vertices).unwrap();
				tmp += " ";
				tmp += &serde_json::to_string(&coef).unwrap();
				tmp += "\n";
				file.write(tmp.to_string().as_bytes()).unwrap();
			}
			counter += 1;
		}
	}

	println!("Complete saving {} representatives!", counter);

///////////////////////////////////////////////////////////////////////////////////////////////////

    Ok(())
}
