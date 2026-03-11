#[derive(Clone, Debug, PartialEq)]
pub struct BigUInt {
    words: Vec<u64>,
}

impl BigUInt {
    pub fn new() -> Self {
        BigUInt { words: Vec::new() }
    }

    pub fn from_u64(n: u64) -> Self {
        if n == 0 {
            BigUInt::new()
        } else {
            BigUInt { words: vec![n] }
        }
    }

    pub fn from_slice(words: &[u64]) -> Self {
        let mut result = BigUInt {
            words: words.to_vec(),
        };
        result.normalize();
        result
    }

    pub fn is_zero(&self) -> bool {
        self.words.is_empty()
    }

    fn normalize(&mut self) {
        while let Some(&last) = self.words.last() {
            if last == 0 {
                self.words.pop();
            } else {
                break;
            }
        }
    }

    pub fn len(&self) -> usize {
        self.words.len()
    }

    pub fn words(&self) -> &[u64] {
        &self.words
    }
}

pub fn next_power_of_two(n: usize) -> usize {
    let mut n = n;
    n -= 1;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    n + 1
}
