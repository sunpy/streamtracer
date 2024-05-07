#[cfg(test)]
mod interp_tests {
    use super::super::interp;
    use numpy::ndarray::array;

    #[test]
    fn test_interp_trilin() {
        let values = array![[[0., 1.], [0., 1.]], [[0., 1.], [0., 1.]]];
        let a = array![0.3, 0.2, 0.4];
        let b = interp::interp_trilinear(&values, &a);
        assert_eq!(b, 0.4);
    }
}
