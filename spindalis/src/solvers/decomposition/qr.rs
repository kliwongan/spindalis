use crate::solvers::SolverError;
use crate::utils::{Arr2D, Arr2DError};

fn norm(matrix: &Arr2D<f64>) -> f64 {
    // Assumes the matrix is 1 dimensional, i.e. a vector
    let mut sum: f64 = 0.0;
    for m in 0..matrix.height {
        for n in 0..matrix.width {
            sum += matrix[m][n].powi(2)
        }
    }

    sum.sqrt()
}

fn submatrix(
    matrix: &Arr2D<f64>,
    hstart: usize,
    hend: usize,
    wstart: usize,
    wend: usize,
) -> Result<Arr2D<f64>, SolverError> {
    // Uses 0-indexing
    if hstart > hend
        || wstart > wend
        || hend >= matrix.height
        || wend >= matrix.width
    {
        return Err(SolverError::XInitOutOfBounds);
    }

    // refactor matrix data
    let mut inner = Vec::<f64>::new();
    for h in hstart..(hend + 1){
        for w in wstart..(wend + 1) {
            inner.push(matrix[h][w]);
        }
    }

    let result = Arr2D::from_flat(inner, 0.0, hend - hstart + 1, wend - wstart + 1);
    if let Ok(result) = result {
        return Ok(result);
    } else {
        println!("Error when creating result");
        return Err(SolverError::XInitOutOfBounds);
    }
}

fn sign(entry: f64) -> f64 {
    if entry < 0.0 {
        return 1.0;
    }
    -1.0
}

pub fn qr_factorization<M>(matrix: M) -> Result<(Arr2D<f64>, Arr2D<f64>), SolverError>
where
    M: TryInto<Arr2D<f64>, Error = Arr2DError>,
{
    // Returns the compressed form of the QR factorization
    let mut matrix: Arr2D<f64> = matrix.try_into()?;
    if matrix.height < matrix.width {
        return Err(SolverError::NonSquareMatrix);
    }

    let mut q: Arr2D<f64> = Arr2D::identity(matrix.height);

    let mut tau: Vec<f64> = Vec::new();
    for j in 0..matrix.width {
        let hvec = submatrix(&matrix, j, matrix.width - 1, j, j).unwrap();
        let normx = norm(&hvec);
        let s = -sign(matrix[j][j]);
        let u1 = matrix[j][j] - s * normx;
        let mut w = hvec / u1;

        // I can't currently think of a better way to mutate the first element
        // of an Arr2D which doesn't have a statically defined size at compile time
        let w_iter = w.rows_mut();
        for row in w_iter.into_iter() {
            for elem in row {
                *elem = 1_f64;
                break;
            }
            break;
        }

        // matrix[j+1:end, j] = w(2:end)
        // let w_iter = w.rows_mut();
        // let mut idx = j + 1;
        // for row in w_iter.into_iter() {
        //     for elem in row {
        //         if idx < matrix.height {
        //             matrix[idx][j] = *elem;
        //             idx += 1;
        //         }
        //     }
        // }

        //matrix[j][j] = s * normx;

        tau.push(-s * u1 / normx);


        // Now modify the original matrix
        for row in j..matrix.height {
            let idx = j + 1;
            let w_iter = w.rows_mut();
            for wrow in w_iter.into_iter() {
                for elem in wrow {
                    if idx < matrix.height {
                        matrix[row][idx] =
                            matrix[row][idx] - (tau[j] * *elem) * (*elem * matrix[row][idx]);
                        q[row][idx] = q[row][idx] - (q[row][idx] * *elem) * (tau[j] * *elem)
                            
                    }
                }
            }
        }
    }

    Ok((q, matrix))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_norm() {
        let vec: Arr2D<f64> = Arr2D::from(&[[-1.0, 3.0, 5.0]]);
        let result = norm(&vec);
        assert_eq!(result, 5.916079783099616);
    }

    #[test]
    fn basic_submatrix() {
        let vec = Arr2D::from(&[[-1.0, 2.0, 3.0], [6.0, 5.0, 4.0]]);
        assert_eq!(vec.shape(), (2, 3));
        let result = submatrix(&vec, 0, 1, 1, 2).unwrap();
        let expected: Arr2D<f64> = Arr2D::from(&[[2.0, 3.0], [5.0, 4.0]]);
        assert_eq!(result, expected);
    }

    #[test]
    fn basic_qr() {
        let matr: Arr2D<f64> = Arr2D::from(&[[-1.0, 3.0], [1.0, 5.0]]);
        let (q, r) = qr_factorization(&matr).unwrap();
        println!("q: {} \n r: {}", q, r);
        //let expected = Arr2D::from(&[[-1.0, 2.0, 3.0], [6.0, 5.0, 4.0]]);
        let expected = q * r;
        //println!("{:?}", result);
        assert_eq!(matr, expected);
    }
}


