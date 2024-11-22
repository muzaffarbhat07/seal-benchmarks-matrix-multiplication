#include "seal/seal.h"
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <seal/ckks.h>
#include <seal/context.h>
#include <seal/encryptor.h>
#include <seal/evaluator.h>
#include <seal/galoiskeys.h>
#include <seal/modulus.h>
#include <vector>
#include <random>
#include <fstream>
#include "seal/util/polyarithsmallmod.h"
#include "matrix_mul.h"
// #include <seal/encryptor.h>
// #include <seal/modulus.h>
// #include <seal/publickey.h>
// #include <seal/batchencoder.h>
// #include <seal/ciphertext.h>
// #include <seal/plaintext.h>

using namespace std;
using namespace seal;
using namespace seal::util;

void print_vector(const vector<double>& vec, size_t size) {
    for (size_t i = 0; i < size; i++) {
        cout << vec[i] << " ";
    }
    cout << endl;
}

void batch_input_operations(Encryptor& encryptor, Decryptor& decryptor, Evaluator& evaluator, CKKSEncoder& encoder, double scale) {
    // Step 3: Encode and encrypt input data
    vector<double> input1 = { 1.0, 2.0, 3.0, 4.0 };
    vector<double> input2 = { 5.0, 6.0, 7.0, 8.0 };
    Plaintext plain1, plain2;
    encoder.encode(input1, scale, plain1);
    encoder.encode(input2, scale, plain2);

    Ciphertext encrypted1, encrypted2;
    encryptor.encrypt(plain1, encrypted1);
    encryptor.encrypt(plain2, encrypted2);

    // Step 4: Ciphertext addition
    Ciphertext encrypted_addition;
    evaluator.add(encrypted1, encrypted2, encrypted_addition);

    // Step 5: Ciphertext scalar multiplication
    double scalar = 3.0;
    Plaintext plain_scalar;
    encoder.encode(vector<double>(input1.size(), scalar), scale, plain_scalar);

    Ciphertext encrypted_scalar_multiplication;
    evaluator.multiply_plain(encrypted1, plain_scalar, encrypted_scalar_multiplication);

    // Step 6: Perform ciphertext multiplication
    Ciphertext encrypted_product;
    evaluator.multiply(encrypted1, encrypted2, encrypted_product);

    // Rescale the result
    evaluator.rescale_to_next(encrypted_product, encrypted_product);

    // Adjust scales and modulus chain if necessary
    encrypted1.scale() = encrypted_product.scale(); // Align scales
    evaluator.mod_switch_to_inplace(encrypted1, encrypted_product.parms_id()); // Align modulus chains

    // Step 7: Decrypt results
    Plaintext plain_result_addition, plain_result_scalar_multiplication, plain_result_ciphertext_multiplication;
    vector<double> result_addition, result_scalar_multiplication, result_ciphertext_multiplication;

    decryptor.decrypt(encrypted_addition, plain_result_addition);
    encoder.decode(plain_result_addition, result_addition);

    decryptor.decrypt(encrypted_scalar_multiplication, plain_result_scalar_multiplication);
    encoder.decode(plain_result_scalar_multiplication, result_scalar_multiplication);

    decryptor.decrypt(encrypted_product, plain_result_ciphertext_multiplication);
    encoder.decode(plain_result_ciphertext_multiplication, result_ciphertext_multiplication);

    // Step 8: Display results
    cout << "Input1: ";
    print_vector(input1, input1.size());

    cout << "Input2: ";
    print_vector(input2, input2.size());

    cout << "Addition Result: ";
    print_vector(result_addition, input1.size());

    cout << "Scalar Multiplication Result: ";
    print_vector(result_scalar_multiplication, input1.size());

    cout << "Ciphertext Multiplication Result: ";
    print_vector(result_ciphertext_multiplication, input1.size());
}

void single_input_operations(Encryptor& encryptor, Decryptor& decryptor, Evaluator& evaluator, CKKSEncoder& encoder, double scale) {
    // Step 3: Encode and encrypt single input numbers
    double num1 = 3.0, num2 = 5.0;

    Plaintext plain1, plain2;
    encoder.encode(num1, scale, plain1); // Encode single number
    encoder.encode(num2, scale, plain2); // Encode single number

    Ciphertext encrypted1, encrypted2;
    encryptor.encrypt(plain1, encrypted1);
    encryptor.encrypt(plain2, encrypted2);

    // Step 4: Perform addition
    Ciphertext encrypted_addition;
    evaluator.add(encrypted1, encrypted2, encrypted_addition);

    // Step 5: Perform scalar multiplication
    double scalar = 4.0;
    Plaintext plain_scalar;
    encoder.encode(scalar, scale, plain_scalar);

    Ciphertext encrypted_scalar_multiplication;
    evaluator.multiply_plain(encrypted1, plain_scalar, encrypted_scalar_multiplication);

    // Step 6: Perform ciphertext multiplication
    Ciphertext encrypted_product;
    evaluator.multiply(encrypted1, encrypted2, encrypted_product);

    // Rescale the result
    evaluator.rescale_to_next(encrypted_product, encrypted_product);

    // Adjust scales and modulus chain if necessary
    encrypted1.scale() = encrypted_product.scale(); // Align scales
    evaluator.mod_switch_to_inplace(encrypted1, encrypted_product.parms_id()); // Align modulus chains

    // Step 6: Decrypt results
    Plaintext plain_result_addition, plain_result_scalar_multiplication, plain_result_ciphertext_multiplication;
    vector<double> result_addition, result_scalar_multiplication, result_ciphertext_multiplication;

    decryptor.decrypt(encrypted_addition, plain_result_addition);
    encoder.decode(plain_result_addition, result_addition); // Decode single number

    decryptor.decrypt(encrypted_scalar_multiplication, plain_result_scalar_multiplication);
    encoder.decode(plain_result_scalar_multiplication, result_scalar_multiplication); // Decode single number

    decryptor.decrypt(encrypted_product, plain_result_ciphertext_multiplication);
    encoder.decode(plain_result_ciphertext_multiplication, result_ciphertext_multiplication);

    // Step 7: Display results
    cout << "Input1: " << num1 << endl;
    cout << "Input2: " << num2 << endl;

    cout << "Addition Result: " << result_addition[0] << endl;
    cout << "Scalar Multiplication Result: " << result_scalar_multiplication[0] << endl;
    cout << "Ciphertext Multiplication Result: " << result_ciphertext_multiplication[0] << endl;
}

void floating_point_operations_using_ckks() {

    // Step 1: Set up encryption parameters and context
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 32768; // N = 2^15;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    int L = 14;
    vector<int> bit_sizes(L, 50);
    bit_sizes[0] = 60;
    // bit_sizes[L - 1] = 53;
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));

    double scale = pow(2.0, 40);
    SEALContext context(parms);

    // Step 2: Create keys and tools
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, secret_key);
    encryptor.set_public_key(public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);

    /*************************BATCH INPUT***********************************/
    cout << "\n/*************************BATCH INPUT********************************/\n";
    batch_input_operations(encryptor, decryptor, evaluator, encoder, scale);

    /*************************SINGLE INPUT***********************************/
    cout << "\n\n/***********************SINGLE INPUT********************************/\n";
    single_input_operations(encryptor, decryptor, evaluator, encoder, scale);
}

void integer_operations_using_bfv() {
    // Step 1: Set up encryption parameters and SEAL context
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 32768; // N = 2^15;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree, sec_level_type::tc128));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 32)); // Plaintext modulus

    SEALContext context(parms);

    // Step 2: Create keys and tools
    KeyGenerator keygen(context);
    PublicKey public_key;
    keygen.create_public_key(public_key);
    SecretKey secret_key = keygen.secret_key();
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    BatchEncoder batch_encoder(context);

    // Step 3: Encode and encrypt integers
    vector<uint64_t> input1(batch_encoder.slot_count(), 3); // Fill all slots with the value 3
    vector<uint64_t> input2(batch_encoder.slot_count(), 5); // Fill all slots with the value 5

    Plaintext plain1, plain2;
    batch_encoder.encode(input1, plain1);
    batch_encoder.encode(input2, plain2);

    Ciphertext encrypted1, encrypted2;
    encryptor.encrypt(plain1, encrypted1);
    encryptor.encrypt(plain2, encrypted2);

    // Step 4: Perform operations
    // Addition
    Ciphertext encrypted_addition;
    evaluator.add(encrypted1, encrypted2, encrypted_addition);

    // Multiplication
    Ciphertext encrypted_multiplication;
    evaluator.multiply(encrypted1, encrypted2, encrypted_multiplication);

    // Scalar multiplication
    Plaintext scalar;
    batch_encoder.encode(vector<uint64_t>(batch_encoder.slot_count(), 4), scalar);
    Ciphertext encrypted_scalar_multiplication;
    evaluator.multiply_plain(encrypted1, scalar, encrypted_scalar_multiplication);

    // Step 5: Decrypt and decode results
    Plaintext plain_result_addition, plain_result_multiplication, plain_result_scalar_multiplication;
    vector<uint64_t> result_addition, result_multiplication, result_scalar_multiplication;

    decryptor.decrypt(encrypted_addition, plain_result_addition);
    batch_encoder.decode(plain_result_addition, result_addition);

    decryptor.decrypt(encrypted_multiplication, plain_result_multiplication);
    batch_encoder.decode(plain_result_multiplication, result_multiplication);

    decryptor.decrypt(encrypted_scalar_multiplication, plain_result_scalar_multiplication);
    batch_encoder.decode(plain_result_scalar_multiplication, result_scalar_multiplication);

    // Step 6: Display results
    cout << "Input1: " << input1[0] << endl;
    cout << "Input2: " << input2[0] << endl;

    cout << "Addition result: " << result_addition[0] << endl;
    cout << "Scalar multiplication result: " << result_scalar_multiplication[0] << endl;
    cout << "Multiplication result: " << result_multiplication[0] << endl;
}

// int main() {
    
//     cout << "\nFloating Point Operations: \n";
//     floating_point_operations_using_ckks();

//     cout << "\n\nInterger Operations: \n\n";
//     integer_operations_using_bfv();

//     return 0;
// }

// Utility function to generate a random matrix
vector<vector<uint64_t>> generateRandomMatrix(int rows, int cols) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> dist(0, pow(2, 32) - 1);

    vector<vector<uint64_t>> matrix(rows, vector<uint64_t>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = dist(gen);
        }
    }
    return matrix;
}

// Encrypt a matrix using BFV
vector<Ciphertext> encryptMatrix(
    const vector<vector<uint64_t>> &matrix, Encryptor &encryptor, BatchEncoder &batchEncoder, Evaluator &evaluator, const PublicKey &public_key, const size_t slot_count) {
    
    vector<Ciphertext> encrypted_matrix;
    Plaintext plain;
    Ciphertext encrypted;

    for (const auto &row : matrix) {
        vector<int64_t> extended_row(row.size(), 0);
        copy(row.begin(), row.end(), extended_row.begin());

        batchEncoder.encode(extended_row, plain);
        encryptor.encrypt(plain, encrypted);
        encrypted_matrix.push_back(encrypted);
    }

    return encrypted_matrix;
}

vector<vector<Ciphertext>> encryptMatrix(
    const vector<vector<uint64_t>> &matrix, Encryptor &encryptor, BatchEncoder &batchEncoder) {
    vector<vector<Ciphertext>> encryptedMatrix;
    for (const auto &row : matrix) {
        vector<Ciphertext> encryptedRow;
        for (const auto &value : row) {
            Plaintext plain;// = Plaintext(to_string(value)); // Convert integer to plaintext
            vector<uint64_t> enc{value};
            batchEncoder.encode(enc, plain);
            Ciphertext encrypted;
            encryptor.encrypt(plain, encrypted);
            encryptedRow.push_back(encrypted);
        }
        encryptedMatrix.push_back(encryptedRow);
    }
    return encryptedMatrix;
}

Ciphertext matrixMultiplyElement(
    const vector<Ciphertext> &rowA, const vector<Ciphertext> &colB, Evaluator &evaluator, RelinKeys &relin_keys) {
    Ciphertext result;
    evaluator.multiply(rowA[0], colB[0], result);
    evaluator.relinearize_inplace(result, relin_keys);
    
    for (size_t k = 1; k < rowA.size(); ++k) {
        Ciphertext temp;
        evaluator.multiply(rowA[k], colB[k], temp);
        evaluator.relinearize_inplace(temp, relin_keys);
        evaluator.add_inplace(result, temp);
    }
    return result;
}

vector<vector<Ciphertext>> matrixMultiply(
    const vector<vector<Ciphertext>> &encryptedA,
    const vector<vector<uint64_t>> &B,
    BatchEncoder &batchEncoder,
    Evaluator &evaluator) {
    size_t rows = encryptedA.size();
    size_t cols = B[0].size();
    size_t common_dim = B.size();

    vector<vector<Ciphertext>> encryptedC(rows, vector<Ciphertext>(cols));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            Ciphertext result;
            for (size_t k = 0; k < common_dim; ++k) {
                Plaintext plain_scalar;
                batchEncoder.encode(vector<uint64_t>{B[k][j]}, plain_scalar);
                Ciphertext mul;
                evaluator.multiply_plain(encryptedA[i][k], plain_scalar, mul);

                if(k == 0) {
                    result = mul;
                } else {
                    evaluator.add_inplace(result, mul);
                }
            }
            encryptedC[i][j] = result;
        }
    }

    return encryptedC;
}

// Perform matrix multiplication
vector<Ciphertext> multiplyMatrices(
    const vector<Ciphertext> &encryptedInputMatrix,
    const vector<vector<uint64_t>> &weightMatrix,
    BatchEncoder &batchEncoder,
    Evaluator &evaluator,
    const RelinKeys &relin_keys,
    const size_t slot_count) {
    
    size_t N = weightMatrix.size();
    size_t R = encryptedInputMatrix.size();
    vector<Ciphertext> result(R);

    // for (size_t i = 0; i < n; ++i) {
    //     Ciphertext sum;
    //     bool initialized = false;

    //     for (size_t j = 0; j < n; ++j) {
    //         // Multiply the encrypted row with the plaintext column weight
    //         vector<int64_t> extended_col(n, weightMatrix[j][i]);
    //         // for (size_t k = 0; k < slot_count; ++k) {
    //         //     extended_col[k] = weightMatrix[j][i];
    //         // }

    //         Plaintext plain_col;
    //         batchEncoder.encode(extended_col, plain_col);

    //         Ciphertext temp;
    //         evaluator.multiply_plain(encryptedInputMatrix[j], plain_col, temp);

    //         // Accumulate the results
    //         if (!initialized) {
    //             sum = temp;
    //             initialized = true;
    //         } else {
    //             evaluator.add_inplace(sum, temp);
    //         }
    //     }

    //     // Relinearize and store the result
    //     evaluator.relinearize_inplace(sum, relin_keys);
    //     result.push_back(sum);
    // }

    for (int i = 0; i < R; ++i) {
        Ciphertext temp_result;
        Plaintext plain_weight;
        batchEncoder.encode(weightMatrix[0], plain_weight);
        evaluator.multiply_plain(encryptedInputMatrix[i], plain_weight, temp_result); // Start with first row of weight
        
        for (int j = 1; j < N; ++j) {
            Ciphertext temp;
            batchEncoder.encode(weightMatrix[j], plain_weight);
            evaluator.multiply_plain(encryptedInputMatrix[i], plain_weight, temp); // Multiply with each row of weight
            evaluator.add_inplace(temp_result, temp); // Add to result
        }
        
        result[i] = temp_result; // Store result for this row
    }

    return result;
}

// Decrypt a matrix
vector<vector<uint64_t>> decryptMatrix(
    const vector<Ciphertext> &encryptedMatrix,
    Decryptor &decryptor,
    BatchEncoder &batchEncoder,
    const size_t slot_count,
    size_t rows) {
    
    vector<vector<uint64_t>> result;
    for (size_t i = 0; i < rows; ++i) {
        Plaintext plain;
        decryptor.decrypt(encryptedMatrix[i], plain);

        vector<int64_t> decoded;
        batchEncoder.decode(plain, decoded);

        vector<uint64_t> row(decoded.begin(), decoded.begin() + 4);
        result.push_back(row);
    }
    return result;
}

vector<vector<uint64_t>> decryptMatrix(
    const vector<vector<Ciphertext>> &encryptedMatrix, Decryptor &decryptor, BatchEncoder &batchEncoder) {
    vector<vector<uint64_t>> decryptedMatrix;
    for (const auto &row : encryptedMatrix) {
        vector<uint64_t> decryptedRow;
        for (const auto &encrypted : row) {
            Plaintext plain;
            decryptor.decrypt(encrypted, plain);

            vector<int64_t> decoded;
            batchEncoder.decode(plain, decoded);

            decryptedRow.push_back(decoded[0]);
            // cout << "decoded size: " << decoded.size() << endl;
        }
        decryptedMatrix.push_back(decryptedRow);
    }
    return decryptedMatrix;
}

vector<int> MM_COEFF_MODULI = {60, 40, 60};
double SCALE = pow(2.0, 40);
int R = 8;
int K = 4;
int C = 4;

void MM_test() {
  long logN = 13;
  size_t poly_modulus_degree = 1 << logN;

  EncryptionParameters parms(scheme_type::ckks);
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, MM_COEFF_MODULI));

  SEALContext context(parms, true, sec_level_type::none);

  KeyGenerator keygen(context);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  keygen.create_public_key(public_key);
  RelinKeys relin_keys;
  keygen.create_relin_keys(relin_keys);
  GaloisKeys galois_keys;

  std::vector<std::uint32_t> rots;
  for (int i = 0; i < logN; i++) {
    rots.push_back((poly_modulus_degree + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
  }
  keygen.create_galois_keys(rots, galois_keys);

  Encryptor encryptor(context, public_key);
  CKKSEncoder encoder(context);
  Evaluator evaluator(context);
  Decryptor decryptor(context, secret_key);

  size_t N = encoder.slot_count() * 2;
  MMEvaluator mme(&encoder, &encryptor, &context, &evaluator, &galois_keys, rots, poly_modulus_degree, N, SCALE);

  std::vector<std::vector<double>> matrix_4096x768 = mme.readMatrix("./data/matrixmul_input_m_128_n_768_k_64_batch_128.txt", 4096, 768);
  std::vector<std::vector<double>> matrix_768x64 = mme.readMatrix("./data/matrix_input_n_768_k_64.txt", 768, 64);

  auto matrix_4096x768_T = mme.transposeMatrix(matrix_4096x768);
  auto matrix_768x64_T = mme.transposeMatrix(matrix_768x64);

  std::vector<std::vector<double>> row_pack;
  std::vector<double> row_ct(poly_modulus_degree, 0.0);
  for (auto i = 0; i < 64 * 768; i++) {
    int row = i / 768;
    int col = i % 768;
    row_ct[i % poly_modulus_degree] = matrix_768x64_T[row][col];
    if (i % poly_modulus_degree == (poly_modulus_degree - 1)) {
      row_pack.push_back(row_ct);
    }
  }

  vector<Ciphertext> res;
  auto start = high_resolution_clock::now();
  mme.matrix_mul(matrix_4096x768_T, row_pack, res);
  auto end = high_resolution_clock::now();
  cout << "[MatMul] 4096x768 x 768x64 takes: " << duration_cast<milliseconds>(end - start).count() << " milliseconds"
       << endl;

//   cout << "Result Matrix: " << endl;
//   for (auto& ct : res) {
//       Plaintext row_pt;
//         vector<double> row;
//         decryptor.decrypt(res[0], row_pt);
//         encoder.decode(row_pt, row);

//         for (auto i = 0; i < 4096; i++) {
//             cout << row[i] << " ";
//         }
//         cout << endl;
//   }

  std::vector<std::vector<double>> matrix_4096x64 = mme.readMatrix("./data/matrix_output_m_128_k_64_batch_128.txt", 4096, 64);
  auto matrix_4096x64_T = mme.transposeMatrix(matrix_4096x64);

  // Calculate the error of the first column
  double average_err = 0.0;
  Plaintext res_pt;
  vector<double> mm_res;
  decryptor.decrypt(res[0], res_pt);
  encoder.decode(res_pt, mm_res);

  for (auto i = 0; i < 4096; i++) {
    average_err += fabs(mm_res[i] / 2.0 - matrix_4096x64_T[0][i]);
  }

  std::cout << "Average error: " << average_err / 4096.0 << std::endl;
}

int main() {
    MM_test();
    return 0;
    // Step 1: Set up SEAL context
    EncryptionParameters params(scheme_type::bfv);
    size_t poly_modulus_degree = 32768;// 2^15       8192;
    params.set_poly_modulus_degree(poly_modulus_degree);
    params.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree, sec_level_type::tc128));
    params.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 32));

    auto context = SEALContext(params);
    KeyGenerator keygen(context);
    PublicKey public_key;
    keygen.create_public_key(public_key);
    auto secret_key = keygen.secret_key();
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    BatchEncoder batchEncoder(context);

    size_t slot_count = batchEncoder.slot_count();

    // Step 2: Generate random matrices
    int N = 4; // Size of the matrices
    auto inputMatrix = generateRandomMatrix(8, N);
    auto weightMatrix = generateRandomMatrix(N, N);

    cout << "Input Matrix (encrypted): " << endl;
    for (const auto &row : inputMatrix) {
        for (const auto &elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }

    cout << "Weight Matrix: " << endl;
    for (const auto &row : weightMatrix) {
        for (const auto &elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }

    // Step 3: Encrypt the input matrix
    cout << "Step 3: Encrypt the input matrix\n";
    auto encryptedInputMatrix = encryptMatrix(inputMatrix, encryptor, batchEncoder);

    // Step 4: Perform matrix multiplication
    cout << "Step 4: Perform matrix multiplication\n";
    auto start_multiplication = chrono::high_resolution_clock::now();
    auto encryptedResult = matrixMultiply(encryptedInputMatrix, weightMatrix, batchEncoder, evaluator);
    auto end_multiplication = chrono::high_resolution_clock::now();
    chrono::duration<double> multiplication_duration = end_multiplication - start_multiplication;

    // Step 5: Decrypt and display the result
    cout << "Step 5: Decrypt and display the result\n";
    auto resultMatrix = decryptMatrix(encryptedResult, decryptor, batchEncoder);

    cout << "Result Matrix: " << endl;
    for (const auto &row : resultMatrix) {
        for (const auto &elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }

    cout << "Time Taken for Multiplication: " << multiplication_duration.count() << " seconds" << endl;
    cout << "slot_count: " << slot_count << endl;
    return 0;
}

