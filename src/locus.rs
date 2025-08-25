use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Both,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
            Strand::Both => write!(f, "."),
        }
    }
}

impl From<&str> for Strand {
    fn from(s: &str) -> Self {
        match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            "." | "both" => Strand::Both,
            _ => Strand::Both, // Default to both for unknown values
        }
    }
}

#[derive(Debug, Clone)]
pub struct Locus {
    pub chromosome: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    pub id: Option<String>,
}

impl Locus {
    pub fn new(
        chromosome: String,
        start: usize,
        end: usize,
        strand: Strand,
        id: Option<String>,
    ) -> Self {
        // Ensure start <= end
        let (start, end) = if start <= end {
            (start, end)
        } else {
            (end, start)
        };

        Self {
            chromosome,
            start,
            end,
            strand,
            id,
        }
    }

    pub fn length(&self) -> usize {
        self.end.saturating_sub(self.start) + 1
    }

    #[allow(dead_code)]
    pub fn overlaps(&self, other: &Locus) -> bool {
        self.chromosome == other.chromosome && self.start <= other.end && other.start <= self.end
    }

    #[allow(dead_code)]
    pub fn contains(&self, other: &Locus) -> bool {
        self.chromosome == other.chromosome && self.start <= other.start && other.end <= self.end
    }

    pub fn extend(&self, upstream: u32, downstream: u32) -> Self {
        let (new_start, new_end) = match self.strand {
            Strand::Minus => (
                self.start.saturating_sub(downstream as usize),
                self.end + upstream as usize,
            ),
            Strand::Plus | Strand::Both => (
                self.start.saturating_sub(upstream as usize),
                self.end + downstream as usize,
            ),
        };

        Self {
            chromosome: self.chromosome.clone(),
            start: if new_start == 0 { 1 } else { new_start },
            end: new_end,
            strand: self.strand,
            id: self.id.clone(),
        }
    }

    #[allow(dead_code)]
    pub fn extend_symmetric(&self, extension: u32) -> Self {
        Self {
            chromosome: self.chromosome.clone(),
            start: self.start.saturating_sub(extension as usize),
            end: self.end + extension as usize,
            strand: self.strand,
            id: self.id.clone(),
        }
    }
}

impl fmt::Display for Locus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}({}):{}:{}",
            self.chromosome, self.strand, self.start, self.end
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_locus_creation() {
        let locus = Locus::new("chr1".to_string(), 100, 200, Strand::Plus, None);
        assert_eq!(locus.start, 100);
        assert_eq!(locus.end, 200);
        assert_eq!(locus.length(), 101);
    }

    #[test]
    fn test_locus_creation_reverse_coords() {
        let locus = Locus::new("chr1".to_string(), 200, 100, Strand::Plus, None);
        assert_eq!(locus.start, 100);
        assert_eq!(locus.end, 200);
    }

    #[test]
    fn test_locus_overlaps() {
        let locus1 = Locus::new("chr1".to_string(), 100, 200, Strand::Plus, None);
        let locus2 = Locus::new("chr1".to_string(), 150, 250, Strand::Plus, None);
        let locus3 = Locus::new("chr1".to_string(), 300, 400, Strand::Plus, None);

        assert!(locus1.overlaps(&locus2));
        assert!(locus2.overlaps(&locus1));
        assert!(!locus1.overlaps(&locus3));
    }

    #[test]
    fn test_locus_extend() {
        let locus = Locus::new("chr1".to_string(), 100, 200, Strand::Plus, None);
        let extended = locus.extend(50, 75);
        assert_eq!(extended.start, 50);
        assert_eq!(extended.end, 275);
    }

    #[test]
    fn test_strand_conversion() {
        assert_eq!(Strand::from("+"), Strand::Plus);
        assert_eq!(Strand::from("-"), Strand::Minus);
        assert_eq!(Strand::from("."), Strand::Both);
        assert_eq!(Strand::from("both"), Strand::Both);
    }
}
