use anyhow::{Context, Result};
use clap::{Arg, Command};
use noodles::{
    bam, sam::alignment::record::cigar::op::Kind as CigarOpKind,
};
use rayon::prelude::*;
use std::{collections::HashMap, path::Path};

mod gff;
mod locus;

use gff::{parse_gff, write_gff, GffRecord};
use locus::{Locus, Strand};

#[derive(Debug, Clone)]
pub struct BamStats {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub uniquely_mapped_reads: Option<u64>,
}

#[derive(Debug, Clone)]
pub struct Config {
    pub bam_file: String,
    pub input_file: String,
    pub output_file: Option<String>,
    pub strand: Strand,
    pub unique_reads: Option<f64>,
    pub density: bool,
    pub floor: u32,
    pub extension: u32,
    pub rpm: bool,
    pub total: bool,
    pub cluster_bin_size: Option<u32>,
    pub matrix_bins: Option<u32>,
    pub include_junction_reads: bool,
}

fn main() -> Result<()> {
    let matches = Command::new("BAM2GFF")
        .version("0.1.0")
        .author("Your Name")
        .about("Maps reads from a BAM file to GFF regions")
        .arg(
            Arg::new("bam")
                .short('b')
                .long("bam")
                .value_name("FILE")
                .help("Sorted BAM file to process")
                .required(true),
        )
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("GFF file to process")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output filename"),
        )
        .arg(
            Arg::new("sense")
                .short('s')
                .long("sense")
                .value_name("STRAND")
                .help("Map to '+', '-', or 'both' strands")
                .default_value("both"),
        )
        .arg(
            Arg::new("unique")
                .short('u')
                .long("unique")
                .value_name("MILLION_READS")
                .help("Number of million uniquely mapping reads for normalization")
                .value_parser(clap::value_parser!(f64)),
        )
        .arg(
            Arg::new("density")
                .short('d')
                .long("density")
                .help("Calculate read density for each region")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("floor")
                .short('f')
                .long("floor")
                .value_name("COUNT")
                .help("Minimum read count threshold for density calculation")
                .default_value("0")
                .value_parser(clap::value_parser!(u32)),
        )
        .arg(
            Arg::new("extension")
                .short('e')
                .long("extension")
                .value_name("BP")
                .help("Extend reads by N base pairs")
                .default_value("200")
                .value_parser(clap::value_parser!(u32)),
        )
        .arg(
            Arg::new("rpm")
                .short('r')
                .long("rpm")
                .help("Normalize density to reads per million")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("total")
                .short('t')
                .long("total")
                .help("Return total read count in region")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("cluster")
                .short('c')
                .long("cluster")
                .value_name("BIN_SIZE")
                .help("Output clustergram with fixed bin size")
                .value_parser(clap::value_parser!(u32)),
        )
        .arg(
            Arg::new("matrix")
                .short('m')
                .long("matrix")
                .value_name("NUM_BINS")
                .help("Output matrix with fixed number of bins")
                .value_parser(clap::value_parser!(u32)),
        )
        .arg(
            Arg::new("junction")
                .short('j')
                .long("junction")
                .help("Include junction reads (reads with 'N' in CIGAR)")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    let config = Config {
        bam_file: matches.get_one::<String>("bam").unwrap().clone(),
        input_file: matches.get_one::<String>("input").unwrap().clone(),
        output_file: matches.get_one::<String>("output").cloned(),
        strand: match matches.get_one::<String>("sense").unwrap().as_str() {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            "both" | "." => Strand::Both,
            _ => return Err(anyhow::anyhow!("Invalid strand specification")),
        },
        unique_reads: matches.get_one::<f64>("unique").copied(),
        density: matches.get_flag("density"),
        floor: *matches.get_one::<u32>("floor").unwrap(),
        extension: *matches.get_one::<u32>("extension").unwrap(),
        rpm: matches.get_flag("rpm"),
        total: matches.get_flag("total"),
        cluster_bin_size: matches.get_one::<u32>("cluster").copied(),
        matrix_bins: matches.get_one::<u32>("matrix").copied(),
        include_junction_reads: matches.get_flag("junction"),
    };

    // Validate inputs
    validate_config(&config)?;

    // Check for BAM index
    check_bam_index(&config.bam_file)?;

    // Parse GFF file
    let gff_records = parse_gff(&config.input_file)?;
    println!("Read {} GFF records", gff_records.len());

    // Get BAM statistics
    // println!("Calculating BAM statistics...");
    // let bam_stats = get_bam_stats(&config.bam_file)?;
    // println!("BAM statistics: {:?}", bam_stats);
    // let mmr = calculate_mmr(&config, &bam_stats)?;
    // println!("Using MMR value of {}", mmr);
    let mmr = 1.0;

    // Process each GFF record
    let results = process_gff_records(&config, &gff_records, mmr)?;

    // Write output
    write_gff(&results, &config)?;

    println!(
        "Output written to: {}",
        config.output_file.as_ref().unwrap_or(&"stdout".into())
    );
    Ok(())
}

fn validate_config(config: &Config) -> Result<()> {
    if !Path::new(&config.bam_file).exists() {
        return Err(anyhow::anyhow!(
            "BAM file does not exist: {}",
            config.bam_file
        ));
    }

    if !Path::new(&config.input_file).exists() {
        return Err(anyhow::anyhow!(
            "Input file does not exist: {}",
            config.input_file
        ));
    }

    if config.cluster_bin_size.is_some() && config.matrix_bins.is_some() {
        return Err(anyhow::anyhow!(
            "Cannot specify both cluster and matrix options"
        ));
    }

    Ok(())
}

fn check_bam_index(bam_file: &str) -> Result<()> {
    let bai_file = format!("{bam_file}.bai");
    if !Path::new(&bai_file).exists() {
        return Err(anyhow::anyhow!(
            "BAM index file not found: {}. Please create an index using 'samtools index'",
            bai_file
        ));
    }
    Ok(())
}

#[allow(dead_code)]
fn get_bam_stats(bam_file: &str) -> Result<BamStats> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_file)
        .with_context(|| format!("Failed to open BAM file: {bam_file}"))?;

    let _header = reader.read_header()?;

    let mut total_reads = 0u64;
    let mut mapped_reads = 0u64;

    for result in reader.records() {
        let record = result?;
        total_reads += 1;

        if !record.flags().is_unmapped() {
            mapped_reads += 1;
        }
    }

    Ok(BamStats {
        total_reads,
        mapped_reads,
        uniquely_mapped_reads: None,
    })
}

#[allow(dead_code)]
fn calculate_mmr(config: &Config, bam_stats: &BamStats) -> Result<f64> {
    if let Some(unique_reads) = config.unique_reads {
        if config.rpm {
            Ok(unique_reads)
        } else {
            Ok(1.0)
        }
    } else if config.rpm {
        Ok(bam_stats.mapped_reads as f64 / 1_000_000.0)
    } else {
        Ok(1.0)
    }
}

fn process_gff_records(
    config: &Config,
    gff_records: &[GffRecord],
    mmr: f64,
) -> Result<Vec<GffRecord>> {
    println!("Processing {} GFF records...", gff_records.len());

    let results: Result<Vec<_>> = gff_records
        .par_iter()
        // .iter()
        .enumerate()
        .map(|(i, record)| {
            if i.is_multiple_of(10000) {
                println!("Processed {i} records");
            }
            process_single_gff_record(config, record, mmr)
        })
        .collect();

    results
}

fn process_single_gff_record(
    config: &Config,
    gff_record: &GffRecord,
    mmr: f64,
) -> Result<GffRecord> {
    let locus = Locus::new(
        gff_record.seqname.clone(),
        gff_record.start,
        gff_record.end,
        gff_record.strand,
        gff_record.attributes.get("ID").cloned(),
    );

    // println!("Locus: {:?}", locus);
    // Create search locus with extension
    let search_locus = locus.extend(config.extension, config.extension);

    // println!("Getting reads");
    // Get reads from BAM
    let reads = get_reads_for_locus(config, &search_locus)?;

    // println!("Extending reads");
    // Extend reads based on strand
    let extended_reads = extend_reads(&reads, config.extension);

    // println!("Separating reads");
    // Separate reads by strand relative to the GFF feature
    let (sense_reads, antisense_reads) = separate_reads_by_strand(&extended_reads, &locus);

    // Calculate result based on requested output type
    let result_value =
        if config.density || config.cluster_bin_size.is_some() || config.matrix_bins.is_some() {
            calculate_density_result(config, &locus, &sense_reads, &antisense_reads, mmr)?
        } else if config.total {
            calculate_total_result(config, &sense_reads, &antisense_reads, mmr)
        } else {
            calculate_read_positions(config, &sense_reads, &antisense_reads)
        };

    // Create new GFF record with result
    let mut new_record = gff_record.clone();
    new_record
        .attributes
        .insert("reads".to_string(), result_value);

    Ok(new_record)
}

fn get_reads_for_locus(config: &Config, locus: &Locus) -> Result<Vec<Locus>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&config.bam_file)
        .with_context(|| format!("Failed to open BAM file: {}", &config.bam_file))?;

    let header = reader.read_header()?;

    // Find the reference sequence ID
    let reference_sequences = header.reference_sequences();
    let chr_id = reference_sequences
        .get_index_of(locus.chromosome.as_bytes())
        .ok_or_else(|| {
            anyhow::anyhow!("Chromosome {} not found in BAM header", locus.chromosome)
        })?;

    let mut reads = Vec::new();

    // println!("Computing region");
    // println!("chr: {}, start: {}, end: {}", locus.chromosome, locus.start, locus.end);
    let region = format!("{}:{}-{}", locus.chromosome, locus.start, locus.end).parse()?;
    // println!("Querying region: {}", region);

    let query = reader.query(&header, &region)?;

    for result in query {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        // Check if read overlaps with our region
        if let Some(ref_seq_id) = record.reference_sequence_id() {
            if ref_seq_id.unwrap() == chr_id {
                if let Some(Ok(alignment_start)) = record.alignment_start() {
                    let start_pos = usize::from(alignment_start);
                    let end_pos = start_pos + record.sequence().len();

                    // Check if read overlaps with our locus
                    if start_pos <= locus.end && end_pos >= locus.start {
                        let read_locus = Locus::new(
                            locus.chromosome.clone(),
                            start_pos,
                            end_pos,
                            if record.flags().is_reverse_complemented() {
                                Strand::Minus
                            } else {
                                Strand::Plus
                            },
                            None,
                        );

                        // Check for junction reads if needed
                        if !config.include_junction_reads && has_junction(&record) {
                            continue;
                        }

                        reads.push(read_locus);
                    }
                }
            }
        }
    }

    Ok(reads)
}

fn has_junction(record: &bam::Record) -> bool {
    record.cigar().iter().any(|op_result| {
        if let Ok(op) = op_result {
            matches!(op.kind(), CigarOpKind::Skip)
        } else {
            false
        }
    })
}

fn extend_reads(reads: &[Locus], extension: u32) -> Vec<Locus> {
    reads
        .iter()
        .map(|read| match read.strand {
            Strand::Plus | Strand::Both => Locus::new(
                read.chromosome.clone(),
                read.start,
                read.end + extension as usize,
                read.strand,
                read.id.clone(),
            ),
            Strand::Minus => Locus::new(
                read.chromosome.clone(),
                read.start.saturating_sub(extension as usize),
                read.end,
                read.strand,
                read.id.clone(),
            ),
        })
        .collect()
}

fn separate_reads_by_strand(reads: &[Locus], gff_locus: &Locus) -> (Vec<Locus>, Vec<Locus>) {
    let mut sense_reads = Vec::new();
    let mut antisense_reads = Vec::new();

    for read in reads {
        match (gff_locus.strand, read.strand) {
            (Strand::Plus, Strand::Plus)
            | (Strand::Plus, Strand::Both)
            | (Strand::Minus, Strand::Minus)
            | (Strand::Minus, Strand::Both)
            | (Strand::Both, _) => {
                sense_reads.push(read.clone());
            }
            (Strand::Plus, Strand::Minus) | (Strand::Minus, Strand::Plus) => {
                antisense_reads.push(read.clone());
            }
        }
    }

    (sense_reads, antisense_reads)
}

fn calculate_density_result(
    config: &Config,
    locus: &Locus,
    sense_reads: &[Locus],
    antisense_reads: &[Locus],
    mmr: f64,
) -> Result<String> {
    let mut coverage = HashMap::new();
    let mut sense_coverage = HashMap::new();
    let mut antisense_coverage = HashMap::new();

    // Fill coverage map
    for read in sense_reads.iter().chain(antisense_reads.iter()) {
        // println!("Read: {:?}", read);
        for pos in read.start..=read.end {
            if read.strand == Strand::Plus || read.strand == Strand::Both {
                *sense_coverage.entry(pos).or_insert(0u32) += 1;
                if pos >= locus.start && pos <= locus.end {
                    *coverage.entry(pos).or_insert(0u32) += 1;
                }
            } else {
                *antisense_coverage.entry(pos).or_insert(0u32) += 1;
                if pos >= locus.start && pos <= locus.end {
                    *coverage.entry(pos).or_insert(0u32) += 1;
                }
            }
        }
    }

    // println!("sense coverage: {:?}", sense_coverage);
    // let mut l: Vec<(&usize, &u32)> = antisense_coverage.iter().collect();
    // l.sort_by_key(|(k,_)| **k);
    // for (k,v) in &l {
    //     println!("{}: {}", k, v);
    // }
    // println!("antisense coverage: {:?}", l);

    // Apply floor filtering
    if config.floor > 0 {
        coverage.retain(|_, &mut count| count > config.floor);
    }

    // Apply coordinate filtering
    coverage.retain(|&pos, _| pos > locus.start && pos < locus.end);

    if let Some(bin_size) = config.cluster_bin_size {
        calculate_clustergram_result(locus, &coverage, bin_size, mmr)
    } else if let Some(num_bins) = config.matrix_bins {
        calculate_matrix_result(locus, &coverage, num_bins, mmr)
    } else {
        // Regular density calculation
        let total_density = coverage.values().sum::<u32>() as f64 / locus.length() as f64;
        let normalized_density = if config.rpm {
            total_density / mmr
        } else {
            total_density
        };
        Ok(format!("{normalized_density:.4}"))
    }
}

fn calculate_clustergram_result(
    locus: &Locus,
    coverage: &HashMap<usize, u32>,
    bin_size: u32,
    mmr: f64,
) -> Result<String> {
    let mut bins = Vec::new();
    let mut current_pos = locus.start;

    while current_pos < locus.end {
        let bin_end = std::cmp::min(current_pos + bin_size as usize, locus.end);
        let bin_coverage: u32 = (current_pos..bin_end)
            .map(|pos| coverage.get(&pos).unwrap_or(&0))
            .sum();

        let bin_density = bin_coverage as f64 / bin_size as f64;
        bins.push(format!("{:.4}", bin_density / mmr));
        current_pos = bin_end;
    }

    Ok(bins.join("\t"))
}

fn calculate_matrix_result(
    locus: &Locus,
    coverage: &HashMap<usize, u32>,
    num_bins: u32,
    mmr: f64,
) -> Result<String> {
    let bin_size = (locus.length() as f64 - 1.0) / num_bins as f64;
    // println!("bin_size: {}", bin_size);
    let mut bins = Vec::new();

    // println!("coverage: {:?}", coverage);

    for i in 0..num_bins {
        let bin_start = locus.start as f64 + (i as f64 * bin_size);
        let mut bin_end = locus.start as f64 + ((i + 1) as f64 * bin_size);
        if bin_end == bin_end.trunc() {
            bin_end -= 0.1; // Adjust to avoid counting a position twice
        }

        // println!("bin_start: {}, bin_end: {}", bin_start, bin_end);

        // (bin_start.ceil() as usize..=bin_end.floor() as usize)
        //     .for_each(|pos| println!("{}: {}", pos, coverage.get(&pos).unwrap_or(&0)));

        // The original Python is not including the bin start position.
        let bin_coverage: u32 = (bin_start.ceil() as usize..=bin_end.floor() as usize)
            .map(|pos| coverage.get(&pos).unwrap_or(&0))
            .sum();

        // println!("bin_coverage: {}", bin_coverage);
        let bin_density = bin_coverage as f64 / bin_size;
        // println!("bin_density: {}", bin_density);
        // Try to match Python output format.
        let mut s = format!("{:.4}", bin_density / mmr);
        s = s.trim_end_matches("0").to_string();

        if (bin_density / mmr) > 0.0 {
            bins.push(s);
        } else {
            bins.push("0.0".to_string());
        }

        // if bin_density == 0.0 {
        //     // To match the original Python implementation, 0 should be rendered with a single digit of precision.
        //     bins.push("0.0".to_string());
        // } else {
        //     bins.push(format!("{:.4}", bin_density / mmr));
        // }
    }
    if locus.strand == Strand::Minus {
        bins.reverse();
    }
    Ok(bins.join("\t"))
}

fn calculate_total_result(
    config: &Config,
    sense_reads: &[Locus],
    antisense_reads: &[Locus],
    mmr: f64,
) -> String {
    let total_count = match config.strand {
        Strand::Plus => sense_reads.len(),
        Strand::Minus => antisense_reads.len(),
        Strand::Both => sense_reads.len() + antisense_reads.len(),
    };

    format!("{:.4}", total_count as f64 / mmr)
}

fn calculate_read_positions(
    config: &Config,
    sense_reads: &[Locus],
    antisense_reads: &[Locus],
) -> String {
    match config.strand {
        Strand::Plus => {
            let positions: Vec<String> = sense_reads.iter().map(|r| r.start.to_string()).collect();
            format!("+:{}", positions.join(","))
        }
        Strand::Minus => {
            let positions: Vec<String> = antisense_reads
                .iter()
                .map(|r| r.start.to_string())
                .collect();
            format!("-:{}", positions.join(","))
        }
        Strand::Both => {
            let sense_positions: Vec<String> =
                sense_reads.iter().map(|r| r.start.to_string()).collect();
            let antisense_positions: Vec<String> = antisense_reads
                .iter()
                .map(|r| r.start.to_string())
                .collect();
            format!(
                "+:{};-:{}",
                sense_positions.join(","),
                antisense_positions.join(",")
            )
        }
    }
}
