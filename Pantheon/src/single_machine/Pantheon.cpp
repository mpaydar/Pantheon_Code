#pragma GCC diagnostic ignored "-Wnarrowing"

#include <cstddef>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>
#include <sstream>
#include <cmath>
#include <ctime>
#include <stack>
#include <unistd.h>
#include <pthread.h>
#include <time.h>
#include<cassert>
#include <fstream> 

#include <iostream>
#include "seal/seal.h"

#include<globals.h>
#include "utils.hpp"
#include <openssl/sha.h>
#include <condition_variable>

using namespace std;
using namespace seal;

#define NUM_COL_THREAD (NUM_COL)
#define NUM_ROW_THREAD NUM_THREAD
#define NUM_PIR_THREAD 4

int TOTAL_MACHINE_THREAD = 6;

#define NUM_EXPANSION_THREAD (TOTAL_MACHINE_THREAD / NUM_COL)
#define NUM_EXPONENT_THREAD (TOTAL_MACHINE_THREAD / (NUM_COL_THREAD * NUM_ROW_THREAD))
// #define NUM_EXPONENT_THREAD (TOTAL_MACHINE_THREAD / (NUM_COL_THREAD ))

typedef  std::vector<seal::Ciphertext> PIRQuery;
typedef  seal::Ciphertext PIRResponse;

struct column_thread_arg {

    int col_id;
    int row_idx;
    Ciphertext *column_result;
    column_thread_arg(int c_id, int r_idx, Ciphertext *res) {
        col_id = c_id;
        row_idx = r_idx;
        column_result = res;
    }
};

struct mult_thread_arg {
    int id;
    int diff;
    Ciphertext *column_result;
    mult_thread_arg(int _id, int _diff, Ciphertext *_column_result) {
        id = _id;
        diff = _diff;
        column_result = _column_result;
    }
};

vector<vector<Plaintext>> db;
mutex mtx; // Mutex for protecting shared state
int stage = 0; // Stage tracker: 0 = expand_query, 1 = process_rows, 2 = process_pir
condition_variable cv; // Condition variable for signaling
size_t slot_count;
size_t row_size;



void sha256(const char *str, int len, unsigned char *dest);
void preallocate_memory();
void *expand_query(void *arg);
void *process_rows(void *arg);
void populate_db();
void *process_columns(void *arg);
void *multiply_columns(void *arg);
void *process_pir(void *arg);
void print_report();
void print_db(const std::vector<std::vector<Plaintext>>& db);

void init_pir_params();
uint32_t get_next_power_of_two(uint32_t number);
uint32_t get_number_of_bits(uint64_t number);
void set_pir_db(std::vector<std::vector<uint64_t> > db);
void pir_encode_db(std::vector<std::vector<uint64_t>> db);
Ciphertext get_sum(Ciphertext *query, uint32_t start, uint32_t end);
vector<uint64_t> rotate_plain(std::vector<uint64_t> original, int index);

int NUM_ROW = 32;
int NUM_COL = 8;
int NUM_THREAD = 1;
vector<Plaintext> masks;
uint64_t number_of_items = 0;

SEALContext *context;

MemoryPoolHandle *column_pools;

Encryptor *encryptor;
Evaluator *evaluator;
BatchEncoder *batch_encoder;
Decryptor *decryptor;
Ciphertext *expanded_query;
Ciphertext *row_result;
Ciphertext *client_query_ct, *one_ct;
Ciphertext *server_query_ct;
Ciphertext *column_results;
Ciphertext *pir_results;

GaloisKeys galois_keys;
RelinKeys relin_keys;

seal::parms_id_type compact_pid;

pthread_t *query_expansion_thread;
int *expansion_thread_id;

pthread_t *row_process_thread;
int *row_thread_id;

pthread_t *pir_thread;
int *pir_thread_id;
//pthread_t *col_process_thread;
int *col_thread_id;

uint64_t *mult_time;
uint64_t *col_processing_time;
uint64_t *row_agg_time;


uint64_t query_gen_time, expansion_time, step1_time, step2_time, total_time;

/****PIR params ******************/

vector<Plaintext> pir_encoded_db;
uint32_t pir_num_obj;
uint32_t pir_obj_size;
uint32_t key_size; //  in bits
uint32_t pir_num_columns_per_obj;
uint32_t pir_num_query_ciphertext;
uint32_t pir_plain_bit_count;
uint32_t pir_db_rows;
/***********************************/

pthread_mutex_t print_lock = PTHREAD_MUTEX_INITIALIZER;

unsigned long *row_thread_processing_time_arr;

int main(int argc, char **argv)
{

    int option;
    const char *optstring = "n:k:s:t:";
    while ((option = getopt(argc, argv, optstring)) != -1)
    {
        switch (option)
        {
        case 'n':
            number_of_items = stoi(optarg);
            break;
        case 's':
            pir_obj_size = stoi(optarg);
            break;
        case 'k':
            key_size = stoi(optarg);
            break;
        case 't':
            NUM_THREAD = stoi(optarg);
            break;

        case '?':
            cout << "error optopt: " << optopt << endl;
            cout << "error opterr: " << opterr << endl;
            return 1;
        }
    }
    if(!number_of_items) {cout<<"Missing -n\n";return 0;}
    if(!NUM_THREAD) {cout<<"Missing -t\n";return 0;}
    if(!key_size) {cout<<"Missing -k\n"; return 0;}
    if(!pir_obj_size) {cout<<"Missing -s\n";return 0;} 

    NUM_COL = (int) ceil(key_size / (2.0 * PLAIN_BIT));

    NUM_ROW = (int) ceil(number_of_items / ((double)(N / 2)));

    init_pir_params();

    chrono::high_resolution_clock::time_point time_start, time_end, total_start, total_end;
    clock_t total_cpu_start, total_cpu_end, cpu_start, cpu_end;
    std::stringstream ss, qss;
    srand(time(NULL));

    mult_time = new uint64_t[NUM_ROW];
    col_processing_time = new uint64_t[NUM_ROW];
    row_agg_time = new uint64_t[NUM_ROW];

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(N);
    parms.set_coeff_modulus(CoeffModulus::Create(N, CT_PRIMES));


    parms.set_plain_modulus(PLAIN_MODULUS);

    context = new SEALContext(parms);
    auto pid = context->first_parms_id();
    uint64_t plain_modulus = parms.plain_modulus().value();
    KeyGenerator keygen(*context);
    SecretKey secret_key = keygen.secret_key();
    keygen.create_relin_keys(relin_keys);


    set<int> rotation_steps;  
    rotation_steps.insert(0);

    for (int i = N / (2 * NUM_COL); i < N / 2; i *= 2) {
        rotation_steps.insert(i);
    }

    for (int i = 1; i < (pir_num_columns_per_obj / 2); i *= 2)
    {
        rotation_steps.insert(-i);
    }

    cout << "Generating Galois keys for steps: ";
    for (int step : rotation_steps) {
        cout << step << " ";
    }
    cout << endl;
    keygen.create_galois_keys(vector<int>(rotation_steps.begin(), rotation_steps.end()), galois_keys);

    encryptor = new  Encryptor(*context, secret_key);
    evaluator = new Evaluator(*context);
    decryptor = new Decryptor(*context,secret_key);

    column_pools = new MemoryPoolHandle[NUM_COL];
    for(int i = 0; i < NUM_COL; i++) {
        column_pools[i] = MemoryPoolHandle::New();
    }
    
    batch_encoder = new BatchEncoder(*context);
    slot_count = batch_encoder->slot_count();
    row_size = slot_count / 2;

    client_query_ct = new Ciphertext(*context);
    server_query_ct = new Ciphertext(*context);
    one_ct = new Ciphertext(*context);

    vector<uint64_t> temp_mat;
    Plaintext temp_pt;

    vector<uint64_t> client_query_mat(slot_count, 0ULL), result_mat, one_mat;
    // Fill the first half with 12
    for (int i = 0; i < slot_count / 2; i++) {
        client_query_mat[i] = 12;
    }

    // Fill the second half with 13
    for (int i = slot_count / 2; i < slot_count; i++) {
        client_query_mat[i] = 13;
    }
    for (int i = 0; i < N ; i++)
    {
        one_mat.push_back(1);
    }
    Plaintext one_pt;
    batch_encoder->encode(one_mat, one_pt);

    encryptor->encrypt_symmetric(one_pt, *one_ct);

    for (int k = 0; k < MOD_SWITCH_COUNT; k++)
    {
        evaluator->mod_switch_to_next_inplace(*one_ct);
    }

    compact_pid = one_ct->parms_id();
    populate_db();

    print_db(db);
    // Debug output to confirm db size and structure
    cout << "DB size after population: " << db.size() << " rows" << endl;
    for (size_t i = 0; i < db.size(); ++i) {
        cout << "Row " << i << " has " << db[i].size() << " columns." << endl;
    }
    vector<vector<uint64_t>> pir_db; //64 bit placeholder for 16 bit plaintext coefficients

    for(int i = 0; i < pir_num_obj; i++) {
        vector<uint64_t> v;
        for(int j = 0; j < (pir_obj_size / 2); j++) {  // 2 bytes each plaintxt slot
            v.push_back(rand() % PLAIN_MODULUS);
        }
        pir_db.push_back(v);
    }

    

    set_pir_db(pir_db);
    cout << "DB population complete!" << endl;


    int desired_index = 0;
    // int val = desired_index + 1;
    // const char str[] = {val & 0xFF, (val >> 8) & 0xFF, (val >> 16) & 0xFF, (val >> 24) & 0xFF, 0};

	// unsigned char hash[SHA256_DIGEST_LENGTH];
    // sha256(str,4, hash);

    // for(int i = 0; i < NUM_COL; i++) {
    //     for(int j = i * (N/(NUM_COL*2)); j < ((i + 1) * (N/(NUM_COL*2))); j++) {
    //         client_query_mat[j] = (uint64_t(hash[4*i]) << 8) + hash[4*i + 1];;
    //         client_query_mat[j + (N/2)] = (uint64_t(hash[4*i + 2]) << 8) + hash[4*i + 3];;
    //     }
    // }
    // client_query_mat[0]=12;

    client_query_mat[desired_index] = 12;
    row_thread_processing_time_arr = new unsigned long[NUM_ROW_THREAD];
    expanded_query = new Ciphertext[NUM_COL];
    pir_results = new Ciphertext[NUM_PIR_THREAD];

    row_result = new Ciphertext[NUM_ROW];

    row_thread_id = new int[NUM_ROW_THREAD];

    for (int i = 0; i < NUM_ROW_THREAD; i++)
    {
        row_thread_id[i] = i;
    }

    col_thread_id = new int[NUM_COL_THREAD];

    for (int i = 0; i < NUM_COL_THREAD; i++)
    {
        col_thread_id[i] = i;
    }

    Plaintext client_query_pt, result_pt;

    time_start = chrono::high_resolution_clock::now();

    batch_encoder->encode(client_query_mat, client_query_pt);
    // cout << "Encoded client_query_pt: " << client_query_pt.to_string() << endl;

    Serializable<Ciphertext> ser_query = encryptor->encrypt_symmetric(client_query_pt);
    // if (!ser_query) {
    //     cerr << "Encryption failed!" << endl;
    // }
    ser_query.save(qss);
    // seal::Ciphertext temp_query; // Create a Ciphertext object
    // temp_query.load(*context, qss); // Deserialize into temp_query
    // Plaintext decrypted_pt;
//    decryptor->decrypt(temp_query, decrypted_pt); // Now decrypt the deserialized ciphertext
    // cout << "Decrypted query: " << decrypted_pt.to_string() << endl;

    time_end = chrono::high_resolution_clock::now();
    query_gen_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();

    printf("query size (Byte): %lu\n", qss.str().size());


    query_expansion_thread = new pthread_t[NUM_COL];
    expansion_thread_id = new int[NUM_COL];

    for (int i = 0; i < NUM_COL; i++)
    {
        expansion_thread_id[i] = i;
    }

    row_process_thread = new pthread_t[NUM_ROW_THREAD];

    pir_thread = new pthread_t[NUM_PIR_THREAD];

    pir_thread_id = new int[NUM_PIR_THREAD];

    for (int i = 0; i < NUM_PIR_THREAD; i++)
    {
        pir_thread_id[i] = i;
    }

    for (int i = 0; i < NUM_COL; i++)
    {
        vector<uint64_t> mat(N, 0ULL);
        Plaintext pt;
        // for (int j = i * (N / (2 * NUM_COL)); j < (i + 1) * (N / (2 * NUM_COL)); j++)
        // {
        //     mat[j] = mat[j + (N / 2)] = 1;
        // }
        for (int j = i * (N / NUM_COL); j < (i + 1) * (N / NUM_COL); j++) 
        {
              mat[j] = 1; // Set the mask for this column
         }
    // Debug: Print the mask content
    // cout << "Mask " << i << " before encoding: ";
    // cout << "Should be in form of [1,1,...,0,0] or [0,0,...,1,1]" << endl; 
    // for (size_t k = 0; k < 4; k++) {
    //     cout << mat[k] << " ";
    // }
    // cout << endl;

        batch_encoder->encode(mat, pt);
        evaluator->transform_to_ntt_inplace(pt, pid);
        masks.push_back(pt);
    }
    total_start = chrono::high_resolution_clock::now();
    total_cpu_start = clock();

    server_query_ct->load(*context, qss);

    total_start = chrono::high_resolution_clock::now();

    time_start = chrono::high_resolution_clock::now();
    cpu_start = clock();
    my_transform_to_ntt_inplace(*context, *server_query_ct, TOTAL_MACHINE_THREAD);
    for (int i = 0; i < NUM_COL; i++)
    {
        if (pthread_create(&(query_expansion_thread[i]), NULL, expand_query, (void *)&(expansion_thread_id[i])))
        {
            printf("Error creating expansion thread");
        }
    }

    for (int i = 0; i < NUM_COL; i++)
    {
        pthread_join(query_expansion_thread[i], NULL);
    }

        // Move to process_rows stage
    {
        lock_guard<mutex> lock(mtx);
        stage = 1;
        cv.notify_all(); // Notify threads waiting for stage 1
    }





    cpu_end = clock();
    time_end = chrono::high_resolution_clock::now();
    expansion_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
    auto query_expansion_cpu_time = double(cpu_end - cpu_start) / double(CLOCKS_PER_SEC);

    time_start = chrono::high_resolution_clock::now();
    cpu_start = clock();

    if(NUM_ROW_THREAD == 1) {
        process_rows((void *)&(row_thread_id[0]));
    }else {
        for (int i = 0; i < NUM_ROW_THREAD; i++)
        {
            if (pthread_create(&(row_process_thread[i]), NULL, process_rows, (void *)&(row_thread_id[i])))
            {
                printf("Error creating processing thread");
            }

        }

        for (int i = 0; i < NUM_ROW_THREAD; i++)
        {
            pthread_join(row_process_thread[i], NULL);
        }
    }

      // Move to process_pir stage
    {
        lock_guard<mutex> lock(mtx);
        stage = 2;
        cv.notify_all(); // Notify threads waiting for stage 2
    }



    time_end = chrono::high_resolution_clock::now();
    step1_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
   
    time_start = chrono::high_resolution_clock::now();

    for (int i = 0; i < NUM_PIR_THREAD; i++)
    {
        if (pthread_create(&(pir_thread[i]), NULL, process_pir, (void *)&(pir_thread_id[i])))
        {
            printf("Error creating PIR processing thread");
        }
    }

    for (int i = 0; i < NUM_PIR_THREAD; i++)
    {
        pthread_join(pir_thread[i], NULL);
    }
    Ciphertext pir_result;
    for(int i = 1; i < NUM_PIR_THREAD; i++) {
        my_add_inplace(*context, pir_results[0], pir_results[i],pir_result);
    }
    cpu_end = clock();
    time_end = chrono::high_resolution_clock::now();
    step2_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
    auto processing_cpu_time = double(cpu_end - cpu_start) / double(CLOCKS_PER_SEC);
    total_time = (chrono::duration_cast<chrono::microseconds>(time_end - total_start)).count();


    Ciphertext final_result = pir_results[0];
    final_result.save(ss);

    total_cpu_end = clock();
    total_end = chrono::high_resolution_clock::now();

    auto latency = (chrono::duration_cast<chrono::microseconds>(total_end - total_start)).count();

    cout << "Result noise budget " << decryptor->invariant_noise_budget(final_result) << endl;
   
    decryptor->decrypt(final_result, result_pt);
    batch_encoder->decode(result_pt, result_mat);

    vector<uint64_t> decoded_response;
    decoded_response = rotate_plain(result_mat, desired_index % row_size);
    bool incorrect_result = false;
    for (int i = 0; i < pir_obj_size / 4 ; i++)
    {
        if ((pir_db[desired_index][i] != decoded_response[i]) || (pir_db[desired_index][i + pir_obj_size / 4] != decoded_response[i + N/2]))
        {
            incorrect_result = true;
            break;
        }
    }

    if (incorrect_result)
    {
        std::cout << "Result is incorrect!" << std::endl<<std::endl;
    }
    else
    {
        std::cout << "Result is correct!" << std::endl<<std::endl;
        print_report();
    }
    return 0;
}

void print_report() {
    
    cout<<"Query expansion time (ms): "<<expansion_time / 1000<<endl;
    cout<<"Equality check time (ms): "<<step1_time / 1000<<endl;
    cout<<"PIR time (ms): "<<step2_time / 1000<<endl;
    cout<<"Total processing time (ms): "<<total_time / 1000<<endl;
}


// vector<string> split(const string &line, char delimiter) {
//     vector<string> tokens;
//     stringstream ss(line);
//     string token;
//     while (getline(ss, token, delimiter)) {
//         tokens.push_back(token);
//     }
//     return tokens;
// }

// // Utility function to check if a string is numeric
// bool is_numeric(const string &str) {
//     return !str.empty() && all_of(str.begin(), str.end(), ::isdigit);
// }





// // Function to populate the database from a .tbl file
// void populate_db() {
//     string filename = "data.tbl"; // Update to the correct path if necessary
//     ifstream file(filename);

//     // Check if the file opened successfully
//     if (!file.is_open()) {
//         cerr << "Error: Could not open the file " << filename << endl;
//         return;
//     }
//     int max_lines = 10;  // Set a reasonable limit for testing
//     int line_number = 0;
//     string line;

//     while (getline(file, line) && line_number < max_lines) {
//         line_number++;
//         cout << "Processing line " << line_number << ": " << line << endl; // Debugging line

//         if (line.empty()) {
//             cout << "Skipping empty line." << endl;
//             continue;
//         }

//         vector<string> row_values = split(line, '|');
//         const size_t expected_columns = 9;

//         // Check for the expected number of columns
//         if (row_values.size() < expected_columns) {
//             cerr << "Warning: Line " << line_number << " does not have the expected number of columns: " << row_values.size() << endl;
//             continue;
//         }

//         vector<Plaintext> row_plaintexts;
//         vector<uint64_t> single_value_vector;

//         // Convert and encode each numeric field after checking with is_numeric
//         if (is_numeric(row_values[0])) {
//             uint64_t item_id = stoull(row_values[0]);
//             single_value_vector = {item_id};
//             Plaintext pt;
//             batch_encoder->encode(single_value_vector, pt);
//             row_plaintexts.push_back(pt);
//         } else {
//             cerr << "Error: Non-numeric item_id on line " << line_number << endl;
//             continue; // Skip this line if item_id is not numeric
//         }

//         // Extract and encode the manufacturer_id
//         size_t pos = row_values[2].find('#');
//         if (pos != string::npos && is_numeric(row_values[2].substr(pos + 1))) {
//             uint64_t manufacturer_id = stoull(row_values[2].substr(pos + 1));
//             single_value_vector = {manufacturer_id};
//             Plaintext pt;
//             batch_encoder->encode(single_value_vector, pt);
//             row_plaintexts.push_back(pt);
//         } else {
//             cerr << "Error: Non-numeric or malformed manufacturer_id on line " << line_number << endl;
//             continue;
//         }

//         // Extract and encode the brand_id
//         pos = row_values[3].find('#');
//         if (pos != string::npos && is_numeric(row_values[3].substr(pos + 1))) {
//             uint64_t brand_id = stoull(row_values[3].substr(pos + 1));
//             single_value_vector = {brand_id};
//             Plaintext pt;
//             batch_encoder->encode(single_value_vector, pt);
//             row_plaintexts.push_back(pt);
//         } else {
//             cerr << "Error: Non-numeric or malformed brand_id on line " << line_number << endl;
//             continue;
//         }

//         // Check and encode the size
//         if (is_numeric(row_values[5])) {
//             uint64_t size = stoull(row_values[5]);
//             single_value_vector = {size};
//             Plaintext pt;
//             batch_encoder->encode(single_value_vector, pt);
//             row_plaintexts.push_back(pt);
//         } else {
//             cerr << "Error: Non-numeric size on line " << line_number << endl;
//             continue;
//         }

//         // Convert and encode the price if valid
//         try {
//             uint64_t price = static_cast<uint64_t>(stod(row_values[7])); // Convert price to integer
//             single_value_vector = {price};
//             Plaintext pt;
//             batch_encoder->encode(single_value_vector, pt);
//             row_plaintexts.push_back(pt);
//         } catch (const invalid_argument& e) {
//             cerr << "Error: Non-numeric price on line " << line_number << endl;
//             continue; // Skip this line if there's an invalid price
//         }

//         // Hash and encode the description field (index 1)
//         unsigned char hash[SHA256_DIGEST_LENGTH];
//         SHA256(reinterpret_cast<const unsigned char*>(row_values[1].c_str()), row_values[1].length(), hash);

//         uint64_t description_hash = 0;
//         for (int i = 0; i < sizeof(uint64_t); ++i) {
//             description_hash = (description_hash << 8) | hash[i];
//         }

//         single_value_vector = {description_hash};
//         Plaintext pt;
//         batch_encoder->encode(single_value_vector, pt);
//         row_plaintexts.push_back(pt);

//         // Add the row to db
//         db.push_back(row_plaintexts);
//     }

//     file.close();

//     // Debugging output for db structure
//     cout << "DB populated with " << db.size() << " rows." << endl;
//     for (size_t i = 0; i < min(db.size(), static_cast<size_t>(5)); i++) {
//         cout << "Row " << i << ": ";
//         for (const auto &pt : db[i]) {
//             cout << pt.to_string() << " ";
//         }
//         cout << endl;
//     }
// }

void print_db(const std::vector<std::vector<Plaintext>>& db)
{
    // Loop over each row in the database
    for (size_t row_idx = 0; row_idx < db.size(); ++row_idx)
    {
        std::cout << "Row " << row_idx << ":" << std::endl;
        
        // Loop over each column in the row
        for (size_t col_idx = 0; col_idx < db[row_idx].size(); ++col_idx)
        {
            // Decode the plaintext to get the values as a vector of uint64_t
            std::vector<uint64_t> decoded_values;
            batch_encoder->decode(db[row_idx][col_idx], decoded_values);
            
            // Print the contents of the decoded vector (first few and last few values)
            std::size_t print_size = 5;
            std::size_t row_size = decoded_values.size() / 2;

            std::cout << "    Column " << col_idx << ": " << std::endl;
            std::cout << "    [";
            for (std::size_t i = 0; i < print_size; i++)
            {
                std::cout << std::setw(3) << std::right << decoded_values[i] << ",";
            }
            std::cout << std::setw(3) << " ...,";

            for (std::size_t i = row_size - print_size; i < row_size; i++)
            {
                std::cout << std::setw(3) << decoded_values[i] << ((i != row_size - 1) ? "," : " ]\n");
            }

            std::cout << "    [";
            for (std::size_t i = row_size; i < row_size + print_size; i++)
            {
                std::cout << std::setw(3) << decoded_values[i] << ",";
            }
            std::cout << std::setw(3) << " ...,";

            for (std::size_t i = 2 * row_size - print_size; i < 2 * row_size; i++)
            {
                std::cout << std::setw(3) << decoded_values[i] << ((i != 2 * row_size - 1) ? "," : " ]\n");
            }

            std::cout << std::endl;
        }
    }
}





// Function to populate the database
void populate_db()
{
    vector<vector<uint64_t> > mat_db;
    for(int i = 0;  i < NUM_ROW * NUM_COL; i++) {
        vector<uint64_t> v(N, 0ULL);
        mat_db.push_back(v);
    }
    unsigned char hash[SHA256_DIGEST_LENGTH];

    for (uint32_t row = 0; row < NUM_ROW * (N/2); row++)
    {
        uint32_t row_in_vector = row % (N/2);
        uint32_t val = row + 1;
        const char str[] = {val & 0xFF, (val >> 8) & 0xFF, (val >> 16) & 0xFF, (val >> 24) & 0xFF, 0};
        sha256(str,4, hash);
        for(int col = 0; col < NUM_COL; col++) {
            int vector_idx = (row / (N/2)) * NUM_COL + col;
            mat_db[vector_idx][row_in_vector] = (uint64_t(hash[4*col]) << 8) + hash[4*col + 1];
            mat_db[vector_idx][row_in_vector + (N/2)] = (uint64_t(hash[4*col + 2]) << 8) + hash[4*col + 3];

        }
    }

    mat_db[0][0]=12;
    for(int i = 0; i < NUM_ROW; i++) {
        vector<Plaintext> row_partition;
        for(int j = 0; j < NUM_COL; j++) {
            Plaintext pt;
            batch_encoder->encode(mat_db[i * NUM_COL + j], pt);
            row_partition.push_back(pt);
        }
        db.push_back(row_partition);
    }
    return;
}



template <typename T>
inline void print_matrix(std::vector<T> matrix, std::size_t row_size)
{
    /*
    We're not going to print every column of the matrix (there are 2048). Instead
    print this many slots from beginning and end of the matrix.
    */
    std::size_t print_size = 5;

    std::cout << std::endl;
    std::cout << "    [";
    for (std::size_t i = 0; i < print_size; i++)
    {
        std::cout << std::setw(3) << std::right << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = row_size - print_size; i < row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ((i != row_size - 1) ? "," : " ]\n");
    }
    std::cout << "    [";
    for (std::size_t i = row_size; i < row_size + print_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = 2 * row_size - print_size; i < 2 * row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ((i != 2 * row_size - 1) ? "," : " ]\n");
    }
    std::cout << std::endl;
}




void *expand_query(void *arg)
{
    vector<uint64_t> rotated_matrix(slot_count, 0ULL);
    vector<uint64_t> Before_rotate_matrix(slot_count, 0ULL);
    vector<uint64_t> expand_matrix(slot_count, 0ULL);

    cout << "Entering Expand Query...." << endl;
    lock_guard<mutex> lock(mtx);
    int id = *((int *)arg);
    Plaintext plain_after_rotate;
    cout << "Entering expand_query function with id = " << id << endl;

    // array of ciphertexts   
    expanded_query[id] = *server_query_ct;

    if (server_query_ct == nullptr) {
        cerr << "server_query_ct is not initialized!" << endl;
    } else {
        cout << "server_query_ct is initialized." << endl;
    }
    my_multiply_plain_ntt(*context, expanded_query[id], masks[id], NUM_EXPANSION_THREAD);
    my_transform_from_ntt_inplace(*context, expanded_query[id], NUM_EXPANSION_THREAD);
    
  
    Ciphertext temp_ct;
    Plaintext before_rotate;
    Plaintext original_plus_rotate;
    Plaintext expandPlain;

    double scale = expanded_query[id].scale();
    temp_ct.scale() = scale;
    
    for (int i = N / (2 * NUM_COL); i < N / 2; i *= 2) {

        temp_ct = expanded_query[id];
        vector<uint64_t> decoded_before_rotate;
        
   
            // Perform rotation
            Ciphertext temp_ct_copy= temp_ct;
            my_rotate_internal(*context, temp_ct, i, galois_keys, column_pools[id], NUM_EXPANSION_THREAD);
            decryptor->decrypt(temp_ct, plain_after_rotate);
            decryptor->decrypt(temp_ct_copy,before_rotate);
            // Before_rotate_matrix
            batch_encoder->decode(plain_after_rotate, rotated_matrix);
            batch_encoder->decode(before_rotate, Before_rotate_matrix);
            cout << "Before Rotation: " << endl;
            print_matrix(Before_rotate_matrix,row_size);

            cout << "After Rotation: " << endl;
            print_matrix(rotated_matrix,row_size);

            Ciphertext expandCipher;
            // my_custom_add_inplace(*context, temp_ct, before_rotate, expandCipher);
            my_add_inplace(*context, temp_ct, temp_ct_copy, expandCipher);

            decryptor->decrypt(expandCipher,expandPlain);
            batch_encoder->decode(expandPlain, expand_matrix);
            expanded_query[id]=expandCipher;

            cout << "After Addition: " << endl;
            print_matrix(expand_matrix,row_size);


        // Print the current state of expanded_query[id]
        cout << "Expanded query after rotation and addition (iteration " << i << "): ";

       
        } 

        cout << endl;
    
    return NULL;
}


std::mutex row_mtx;

void *process_rows(void *arg)
{
    lock_guard<mutex> lock(mtx);

    int id = *((int *)arg);
    cout << "Entering Process Rows function with id = " << id << endl;
    // Initialize column_results and pir_results with SEALContext
    vector<Ciphertext> column_results(NUM_COL, Ciphertext(*context));
    vector<Ciphertext> pir_results(NUM_PIR_THREAD, Ciphertext(*context));

    vector<column_thread_arg> column_args;
    vector<mult_thread_arg> mult_args;

    for (int i = 0; i < NUM_COL_THREAD; i++) {
        column_args.push_back(column_thread_arg(i, id, column_results.data()));
    }
    for (int i = 0; i < NUM_COL; i++) {
        mult_args.push_back(mult_thread_arg(i, 1, column_results.data()));
    }
    int num_row_per_thread = NUM_ROW / NUM_ROW_THREAD;
    int start_idx = num_row_per_thread * id;
    int end_idx = start_idx + num_row_per_thread;

    pthread_t col_process_thread[NUM_COL_THREAD];
    pthread_t col_mult_thread[NUM_COL_THREAD];

    for (int row_idx = start_idx; row_idx < end_idx; row_idx++) {
        for (int i = 0; i < NUM_COL_THREAD; i++) {
            column_args[i].row_idx = row_idx;

            if (pthread_create(&(col_process_thread[i]), NULL, process_columns, (void *)&(column_args[i]))) {
                printf("Error creating column processing thread");
            }
        }

        for (int i = 0; i < NUM_COL_THREAD; i++) {
            pthread_join(col_process_thread[i], NULL);
        }

        for (int diff = 2; diff <= NUM_COL; diff *= 2) {
            for (int i = 0; i < mult_args.size(); i++) {
                mult_args[i].diff = diff;
            }
            for (int i = 0; i < NUM_COL; i += diff) {
                if (pthread_create(&(col_mult_thread[i]), NULL, multiply_columns, (void *)&(mult_args[i]))) {
                    printf("Error creating column processing thread");
                }
            }

            for (int i = 0; i < NUM_COL; i += diff) {
                pthread_join(col_mult_thread[i], NULL);
            }
        }

        Ciphertext temp_ct = column_results[0];
        my_conjugate_internal(*context, temp_ct, galois_keys, column_pools[0], TOTAL_MACHINE_THREAD / NUM_ROW_THREAD);
        bfv_multiply(*context,*evaluator ,column_results[0], temp_ct,relin_keys);
        my_relinearize_internal(*context, column_results[0], relin_keys, 2, MemoryManager::GetPool(), TOTAL_MACHINE_THREAD / NUM_ROW_THREAD);
        my_transform_to_ntt_inplace(*context, column_results[0], TOTAL_MACHINE_THREAD);
        row_result[row_idx] = column_results[0];
    // Plaintext row_plain;
    // vector<uint64_t> row_results(slot_count, 0ULL);
    // decryptor->decrypt(row_result[row_idx],row_plain);
    // batch_encoder->decode(row_plain,row_results);
    // print_matrix(row_results,row_size);
    }
   




    cout << "Successfully processed rows..." << endl;
    return nullptr;
}


void *process_pir(void *arg) {
    lock_guard<mutex> lock(mtx);
    cout << "Entering process pir..." << endl;
    int my_id = *((int *)arg);
    int column_per_thread = (pir_num_columns_per_obj / 2) / NUM_PIR_THREAD;
    int start_idx = my_id * column_per_thread;
    int end_idx = start_idx + column_per_thread - 1;
    pir_results[my_id] = get_sum(row_result, start_idx, end_idx);

    int mask = 1;
    while(mask <= start_idx) {
        if(start_idx & mask) {
            my_rotate_internal(*context, pir_results[my_id] , -mask, galois_keys, MemoryManager::GetPool(), TOTAL_MACHINE_THREAD/NUM_PIR_THREAD);
        }
        mask <<= 1;

    }

    cout << "successfully done " << endl;
    return nullptr;
}


void *multiply_columns(void *arg) {
    cout << "Entering Multiply Columns Function..." << endl;
    vector<uint64_t> column_aggregation_result(slot_count, 0ULL);
    vector<uint64_t> column_before_result(slot_count, 0ULL);

    mult_thread_arg mult_arg = *((mult_thread_arg *)arg);
    Ciphertext *column_results = mult_arg.column_result;
    int id = mult_arg.id;
    int diff = mult_arg.diff;
    int num_threads = TOTAL_MACHINE_THREAD / (NUM_COL / diff);
    Ciphertext aggregateColumnResult ;
    cout << "Columns Before Multipication: " << endl;
    Plaintext columnBeforeMult;
    decryptor->decrypt(column_results[id],columnBeforeMult);
    batch_encoder->decode(columnBeforeMult,column_before_result);
    print_matrix(column_before_result,row_size);

    bfv_multiply(*context,*evaluator ,column_results[id], column_results[id + (diff/2)],relin_keys);

    // my_bfv_multiply(*context, column_results[id], column_results[id + (diff/2)], column_pools[id], num_threads);
    my_relinearize_internal(*context, column_results[id] ,relin_keys, 2, column_pools[id], num_threads); 
    cout << "Column After Multipication: "<< endl;
    Plaintext col_result;
    decryptor->decrypt(column_results[id],col_result);
    batch_encoder->decode(col_result,column_aggregation_result);
    print_matrix(column_aggregation_result,row_size);
    cout << "Multiplying was done successfully..." << endl;


    return nullptr;
}



void *process_columns(void *arg) {
    vector<uint64_t> query_subtract_dbValue_results(slot_count, 0ULL);
    vector<uint64_t> incoming_quries(slot_count, 0ULL);
    vector<uint64_t> incoming_db_values(slot_count, 0ULL);
    vector<uint64_t> outgoing_results(slot_count, 0ULL);

    std::vector<seal::Ciphertext> subtraction_results(db.size());
    vector<uint64_t> v_integers(db.size(), 0ULL);

    column_thread_arg col_arg = *((column_thread_arg *)arg);
    vector<unsigned long> exp_time;
    unsigned long exponent_time = 0;
    int num_col_per_thread = NUM_COL / NUM_COL_THREAD;
    int start_idx = num_col_per_thread * col_arg.col_id;
    int end_idx = start_idx + num_col_per_thread;
    for(int i = start_idx; i < end_idx; i++) {
        Ciphertext sub;
        Ciphertext prod;

        evaluator->sub_plain(expanded_query[i], db[col_arg.row_idx][0], sub);
        Plaintext incoming_expand_query;
        decryptor->decrypt(expanded_query[i],incoming_expand_query);
        batch_encoder->decode(incoming_expand_query,incoming_quries);
        cout << "Incoming Expanded Quries" << endl;
        print_matrix(incoming_quries,row_size);

        batch_encoder->decode(db[col_arg.row_idx][i],incoming_db_values);
        cout << "Incoming db values: " << endl;
        print_matrix(incoming_db_values,row_size);
      



        Plaintext exapnd_plain;
        decryptor->decrypt(sub,exapnd_plain);
        batch_encoder->decode(exapnd_plain,query_subtract_dbValue_results);
        cout << "Subtraction Results (slot " << i << "):" << endl;
        print_matrix(query_subtract_dbValue_results,row_size);

        // subtraction_results[i]=query_subtract_dbValue_results[0];
        // print_matrix(query_subtract_dbValue_results,row_size);
        for (int k = 0; k < 16; k++){
            my_bfv_square(*context, sub, column_pools[i], NUM_EXPONENT_THREAD);
            my_relinearize_internal(*context, sub, relin_keys, 2, column_pools[i], NUM_EXPONENT_THREAD);
        }
        for(int k = 0; k < MOD_SWITCH_COUNT; k++) {
            my_mod_switch_scale_to_next(*context, sub, sub, column_pools[i], NUM_EXPONENT_THREAD);
            cout << "MOD: " << k << endl;
            Plaintext one_subtract_sub;
            decryptor->decrypt(sub,one_subtract_sub);
            batch_encoder->decode(one_subtract_sub,outgoing_results);
            cout << "Generated Results from my_bfc_square: " << endl;
            print_matrix(outgoing_results,row_size);
        }
        evaluator->sub(*one_ct, sub, (col_arg.column_result)[i]);




    }
   
    cout << endl;
    cout << "Processed Columns Successsfully" << endl;
    return nullptr;
}

/*****************************PIR functions ************************************/
/******************************************************************************/

Ciphertext get_sum(Ciphertext *query, uint32_t start, uint32_t end)
{
    vector<uint64_t> pir_column_result(slot_count, 0ULL);
    vector<uint64_t> query_result(slot_count, 0ULL);

    Plaintext client_query;
    // my_transform_from_ntt_inplace(*context, *query, TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);
    // decryptor->decrypt(*query,client_query);
    // batch_encoder->decode(client_query,query_result);
    // cout << "Client Query: " << endl;
    // print_matrix(query_result,row_size);

    seal::Ciphertext result;
    Ciphertext sum_result;
    if (start != end)
    {
        int count = (end - start) + 1;
        int next_power_of_two = get_next_power_of_two(count);
        int mid = next_power_of_two / 2;
        seal::Ciphertext left_sum = get_sum(query, start, start + mid - 1);
        seal::Ciphertext right_sum = get_sum(query, start + mid, end);
        my_rotate_internal(*context, right_sum, -mid, galois_keys, column_pools[0], TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);
        my_add_inplace(*context, left_sum, right_sum,sum_result);
        return sum_result;
    }
    else
    {

        Plaintext decrypted_pir;
        Ciphertext column_sum = query[0];
      
        // my_transform_from_ntt_inplace(*context, column_sum, TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);
       

        Ciphertext temp_ct;
        Ciphertext result_addition;
        my_multiply_plain_ntt(*context, column_sum, pir_encoded_db[pir_num_query_ciphertext * start], TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);

        for (int j = 1; j < pir_num_query_ciphertext; j++)
        {
            temp_ct = query[j];
            my_multiply_plain_ntt(*context, temp_ct, pir_encoded_db[pir_num_query_ciphertext * start + j], TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);
            my_add_inplace(*context, column_sum, temp_ct, result_addition);
            // decryptor->decrypt(column_sum, decrypted_pir);
            // batch_encoder->decode(decrypted_pir, pir_column_result);
            // cout << "PIR column_sum: " << j<< endl;
            // print_matrix(pir_column_result,row_size);

        }
        
     
            my_transform_from_ntt_inplace(*context, column_sum, TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);


        return column_sum;
    }
}

uint32_t get_next_power_of_two(uint32_t number)
{
    if (!(number & (number - 1)))
    {
        return number;
    }

    uint32_t number_of_bits = get_number_of_bits(number);
    return (1 << number_of_bits);
}


uint32_t get_number_of_bits(uint64_t number)
{
    uint32_t count = 0;
    while (number)
    {
        count++;
        number /= 2;
    }
    return count;
}


void set_pir_db(std::vector<std::vector<uint64_t> > db)
{
    assert(db.size() == pir_num_obj);
    std::vector<std::vector<uint64_t> > extended_db(pir_db_rows);
    for(int i = 0; i < pir_db_rows; i++) {
        extended_db[i] = std::vector<uint64_t>(N, 1ULL);
    }
    int row_size = N/2;

    for(int i = 0; i < pir_num_obj;i++) {
        std::vector<uint64_t> temp = db[i];

        int row = (i / row_size);
        int col = (i % row_size);
        for (int j = 0; j < pir_num_columns_per_obj / 2; j++)
        {
            extended_db[row][col] = temp[j];
            extended_db[row][col+row_size] = temp[j+(pir_num_columns_per_obj / 2)];
            row += pir_num_query_ciphertext;
        }

    }   
    pir_encode_db(extended_db);
    return;
}

void pir_encode_db(std::vector<std::vector<uint64_t>> db)
{
    pir_encoded_db = std::vector<seal::Plaintext>(db.size());
    for (int i = 0; i < db.size(); i++)
    {
        batch_encoder->encode(db[i], pir_encoded_db[i]);
        evaluator->transform_to_ntt_inplace(pir_encoded_db[i], compact_pid);

    }
}

vector<uint64_t> rotate_plain(std::vector<uint64_t> original, int index)
{
    int sz = original.size();
    int row_count = sz / 2;
    std::vector<uint64_t> result(sz);
    for (int i = 0; i < row_count; i++)
    {
        result[i] = original[(index + i) % row_count];
        result[row_count + i] = original[row_count + ((index + i) % row_count)];
    }

    return result;
}

void init_pir_params()
{
    pir_plain_bit_count = 16;
    pir_num_obj = ((N/2) * NUM_ROW);
    pir_num_query_ciphertext = ceil(pir_num_obj / (double)(N/2));
    pir_num_columns_per_obj = 2 * (ceil(((pir_obj_size/2) * 8) / (float)(pir_plain_bit_count)));
    pir_db_rows = ceil(pir_num_obj / (double)N) * pir_num_columns_per_obj;

    return;

}

void sha256(const char *str, int len, unsigned char *dest)
{    
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, str, len);
    SHA256_Final(dest, &sha256);
}