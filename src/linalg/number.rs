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
        let mut mult: i64 = 100;
        let mut base = self.clone() * (mult as f64);

        while base.fract().abs() > f64::EPSILON {
            base -= self;
            mult -= 1;
        }

        Fraction {
            numerator: base / base,
            denominator: mult as f64 / base,
        }
    }
}
