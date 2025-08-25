use anyhow::{Context, Result};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use crate::{locus::Strand, Config};

#[derive(Debug, Clone)]
pub struct GffRecord {
    pub seqname: String,
    pub source: String,
    #[allow(dead_code)]
    pub feature: String,
    pub start: usize,
    pub end: usize,
    pub score: Option<f64>,
    pub strand: Strand,
    pub frame: Option<u8>,
    pub attributes: HashMap<String, String>,
}

impl GffRecord {
    #[allow(clippy::too_many_arguments, dead_code)]
    pub fn new(
        seqname: String,
        source: String,
        feature: String,
        start: usize,
        end: usize,
        score: Option<f64>,
        strand: Strand,
        frame: Option<u8>,
        attributes: HashMap<String, String>,
    ) -> Self {
        Self {
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            frame,
            attributes,
        }
    }

    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.trim().split('\t').collect();

        if fields.len() < 8 {
            return Err(anyhow::anyhow!("Invalid GFF line: insufficient fields"));
        }

        let seqname = fields[0].to_string();
        let source = fields[1].to_string();
        let feature = fields[2].to_string();

        let start: usize = fields[3]
            .parse()
            .with_context(|| format!("Invalid start coordinate: {}", fields[3]))?;

        let end: usize = fields[4]
            .parse()
            .with_context(|| format!("Invalid end coordinate: {}", fields[4]))?;

        let score = if fields[5] == "." {
            None
        } else {
            Some(
                fields[5]
                    .parse()
                    .with_context(|| format!("Invalid score: {}", fields[5]))?,
            )
        };

        let strand = Strand::from(fields[6]);

        let frame = if fields[7] == "." {
            None
        } else {
            Some(
                fields[7]
                    .parse()
                    .with_context(|| format!("Invalid frame: {}", fields[7]))?,
            )
        };

        let mut attributes = HashMap::new();
        if fields.len() > 8 && !fields[8].is_empty() && fields[8] != "." {
            // Parse GFF3-style attributes (key=value;key=value)
            if fields[8].contains('=') {
                for attr in fields[8].split(';') {
                    if let Some((key, value)) = attr.split_once('=') {
                        attributes.insert(key.trim().to_string(), value.trim().to_string());
                    }
                }
            } else {
                // Simple case: treat the whole field as an ID
                attributes.insert("ID".to_string(), fields[8].to_string());
            }
        }

        // Set a default ID if none exists
        if !attributes.contains_key("ID") {
            let id = format!("{seqname}:{start}-{end}");
            attributes.insert("ID".to_string(), id);
        }

        Ok(Self {
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            frame,
            attributes,
        })
    }

    pub fn to_gff_line(&self) -> String {
        let _score_str = self
            .score
            .map(|s| s.to_string())
            .unwrap_or_else(|| ".".to_string());

        let _frame_str = self
            .frame
            .map(|f| f.to_string())
            .unwrap_or_else(|| ".".to_string());

        let attributes_str = if self.attributes.is_empty() {
            ".".to_string()
        } else {
            // self.attributes
            //     .iter()
            //     .map(|(k, v)| format!("{}={}", k, v))
            //     .collect::<Vec<_>>()
            //     .join(";")
            self.attributes
                .get("reads")
                .cloned()
                .unwrap_or_else(|| ".".to_string())
        };

        format!(
            "{}\t{}({}):{}-{}\t{}",
            self.source, self.seqname, self.strand, self.start, self.end, attributes_str
        )
    }
}

pub fn parse_gff<P: AsRef<Path>>(path: P) -> Result<Vec<GffRecord>> {
    let file = File::open(&path)
        .with_context(|| format!("Failed to open GFF file: {}", path.as_ref().display()))?;

    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line {}", line_num + 1))?;

        // Skip comment lines and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let record = GffRecord::from_line(&line)
            .with_context(|| format!("Failed to parse line {}: {}", line_num + 1, line))?;

        records.push(record);
    }

    Ok(records)
}

pub fn write_gff(records: &[GffRecord], config: &Config) -> Result<()> {
    let file = config
        .output_file
        .clone()
        .unwrap_or_else(|| "matrix.txt".into());
    let path = Path::new(&file);
    let file = File::create(path)
        .with_context(|| format!("Failed to create output file: {}", &path.display()))?;

    let mut writer = BufWriter::new(file);

    // Write GFF header
    // writeln!(writer, "##gff-version 3")?;

    let mut v = vec!["GENE_ID".to_string(), "locusLine".to_string()];
    // TODO: pass the bin count and the bam name in to generate the header
    let p = Path::new(&config.bam_file);
    for i in 1..=config.matrix_bins.unwrap_or(50) {
        v.push(format!(
            "bin_{}_{}",
            i,
            p.file_name().unwrap_or_default().to_string_lossy()
        ));
    }
    writeln!(writer, "{}", v.join("\t"))?;

    for record in records {
        writeln!(writer, "{}", record.to_gff_line())?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gff_record_from_line() {
        let line = "chr1\tsource\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=test_gene";
        let record = GffRecord::from_line(line).unwrap();

        assert_eq!(record.seqname, "chr1");
        assert_eq!(record.start, 1000);
        assert_eq!(record.end, 2000);
        assert_eq!(record.strand, Strand::Plus);
        assert_eq!(record.attributes.get("ID"), Some(&"gene1".to_string()));
    }

    #[test]
    fn test_gff_record_to_line() {
        let mut attributes = HashMap::new();
        attributes.insert("ID".to_string(), "gene1".to_string());

        let record = GffRecord::new(
            "chr1".to_string(),
            "source".to_string(),
            "gene".to_string(),
            1000,
            2000,
            None,
            Strand::Plus,
            None,
            attributes,
        );

        let line = record.to_gff_line();
        assert!(line.contains("chr1\tsource\tgene\t1000\t2000\t.\t+\t.\tID=gene1"));
    }
}
