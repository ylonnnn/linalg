use std::fmt;

#[derive(Debug)]
pub struct Fraction {
    numerator: i64,
    denominator: i64,
}

impl fmt::Display for Fraction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}

pub trait Real {
    fn gcd(&self, other: &Self) -> Option<i64>;
    fn fraction(&self) -> Fraction;
}

impl Real for f64 {
    fn gcd(&self, other: &Self) -> Option<i64> {
        if *self == 0_f64 || *other == 0_f64 {
            let max = self.max(other.clone()) as i64;
            if max == 0 { None } else { Some(max.abs()) }
        } else {
            let mut a = self.abs() as i64;
            let mut b = other.abs() as i64;

            while b != 0 {
                let r = a % b;
                a = b;
                b = r;
            }

            Some(a)
        }
    }

    fn fraction(&self) -> Fraction {
        let mut mult: i64 = 100;
        let mut base = self.clone() * (mult as f64);

        while base.fract().abs() > f64::EPSILON {
            base -= self;
            mult -= 1;
        }

        let nume = base;
        let denom = mult as f64;

        let gcd = nume.gcd(&denom).unwrap_or(1);

        Fraction {
            numerator: nume as i64 / gcd,
            denominator: denom as i64 / gcd,
        }
    }
}
