use std::fmt;

#[derive(Debug)]
pub struct Fraction {
    numerator: f64,
    denominator: f64,
}

impl fmt::Display for Fraction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.numerator, self.denominator)
    }
}

pub trait Real {
    fn fraction(&self) -> Fraction;
}

impl Real for f64 {
    fn fraction(&self) -> Fraction {
        todo!()
    }
}
