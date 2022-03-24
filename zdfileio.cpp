#include "zdfileio.hpp"

int get_file_lines(std::fstream &fs)
{
    // Count number of lines from fs (not necessarily start from std::ios_base::beg)

    int pos = fs.tellg();
    int n = 0;
    char str[100000];
    while (fs.getline(str, 100000))
    {
        if (strlen(str) != 0)
            n = n + 1;
    }

    fs.seekg(pos, std::ios_base::beg);
    // Back to starting point.

    return n;
}

int get_file_cols(std::fstream &fs)
{
    int pos = fs.tellg();
    int n = 0;
    char ele;
    bool is_space = true;
    while (fs.get(ele))
    {
        if (ele == ' ' || ele == '\t')
            is_space = true;
        else if (ele == '\n')
            break;
        else
        {
            if (is_space)
                n++;
            is_space = false;
        }
    }

    fs.seekg(pos, std::ios_base::beg);
    return n;
}

String make_stdname(String s, std::map<String, String> &paras)
{
    String Ans = paras["-path"];
    Ans += s;
    if (paras["-post"] != String("none"))
    {
        Ans += "_";
        Ans += paras["-post"];
    }
    Ans += ".out";
    return Ans;
}

String get_path(String fname)
{
    String path = fname;
    std::string temp = (char *)path;
    size_t loc_slash = temp.find_last_of('/');
    if (loc_slash != string::npos)
        path = path(0, loc_slash + 1);
    else
        path = "";
    return path;
}

bool read_subset(String fname, Array<size_t> &subset)
{
    std::ifstream fin;
    fin.open(fname);
    if (!fin.is_open())
        return false;
    size_t n = 0;
    fin >> n;

    size_t temp = 0;
    subset.resize(n, 0);
    subset.erase();
    for (int i = 0; i < n; i++)
    {
        fin >> temp;
        subset.push(temp);
    }
    return true;
}

void read_paras_from_csv(String fname, map<String, String> &paras, bool allow_insert)
{
    char wmsg[200];
    char ele;

    char line[200];
    int pos = 0;

    String key, value;
    bool is_value = false;

    std::string split = " ;,\n\t\r";
    std::string comment = "#%";

    std::ifstream fin;
    fin.open(fname, std::ios::in);
    if (!fin)
    {
        sprintf(wmsg, "fail to open %s file", (char *)fname);
        std::cout << "Error! file : zdfileio.hpp / read_paras_from_csv(String fname, map<String, String> &paras)\n"
                  << wmsg
                  << ".\n";
    }

    // Get the first line.
    while (fin.get(ele))
        if (ele == '\n')
            break;

    while (fin.get(ele))
    {
        if (comment.find(ele) < 2) // % or # enoucnted, the rest of this line get commented out.
        {
            while (fin.get(ele))
                if (ele == '\n')
                    break;
            continue;
        }

        if (split.find(ele) == std::string::npos)
        {
            if (ele != '"')
                line[pos++] = ele;
        }
        else if (pos != 0)
        {
            line[pos] = '\0';
            if (is_value)
            {
                value = String(line);
                auto it = paras.find(key);
                if (it != paras.end())
                {
                    if (it->second != value)
                    {
                        sprintf(wmsg, "Warning! Changing default parameters when reading %s file:\t ", (char *)fname);
                        std::cout << wmsg << key << "\t" << it->second << " ==> " << value << "\n";
                        it->second = value;
                    }
                }
                else
                {
                    if (allow_insert)
                        paras[key] = value;
                    else
                    {
                        sprintf(wmsg, "Error! Unsupported keyword encounted when reading %s file:\t ", (char *)fname);
                        std::cout << wmsg << key << ".\n"
                                  << "Keyword discarded!\n";
                    }
                }

                is_value = false;
            }
            else
            {
                key = String(line);
                is_value = true;
            }
            pos = 0;
            line[pos] = '\0';
        }
    }

    if (paras.find((String) "-f") != paras.end())
        paras["-path"] = get_path(paras[String("-f")]);
}

map<String, String> read_paras(int argc, char *argv[], int key_size, String *default_paras, String *options)
{
    std::cout << "loading all parameters\n";
    for (int i = 1; i < argc; i++)
    {
        for (int j = 0; j < key_size; j++)
        {
            if ((String)argv[i] == options[j] && i + 1 < argc && argv[i + 1][0] != '-')
            {
                i++;
                default_paras[j] = argv[i];
                std::cout << default_paras[j] << "\t" << argv[i] << "\n";
                break;
            }
        }
    }
    std::cout << "\n";
    map<String, String> paras;
    for (int i = 0; i < key_size; i++)
    {
        paras[options[i]] = default_paras[i];
        std::cout << options[i] << "\t" << default_paras[i] << "\n";
    }

    paras["-path"] = get_path(paras[String("-f")]);
    std::cout << "Done!\n";
    return paras;
}

void binary_print_smat(std::ostream &output, SparseMatrix *sb2t_mat)
{
    Array<PRECISION> *val_c_ptr;
    Array<size_t> *row_ind_c_ptr;
    PRECISION *val_ptr;
    size_t *row_ind_ptr;
    size_t container_size;

    size_t r, c;
    r = sb2t_mat->get_row();
    c = sb2t_mat->get_col();

    output.write((char *)&r, sizeof(size_t));
    output.write((char *)&c, sizeof(size_t));

    size_t col = sb2t_mat->get_col();

    for (size_t i = 0; i < col; i++)
    {
        row_ind_c_ptr = sb2t_mat->get_CCS_ind_c_ptr(i);
        val_c_ptr = sb2t_mat->get_CCS_val_c_ptr(i);

        row_ind_ptr = row_ind_c_ptr->get_vec();
        val_ptr = val_c_ptr->get_vec();

        container_size = row_ind_c_ptr->get_size();

        output.write((char *)&container_size, sizeof(size_t));
        size_t ind;
        PRECISION val;
        for (size_t j = 0; j < container_size; j++)
        {
            ind = row_ind_ptr[j];
            val = val_ptr[j];
            output.write((char *)&ind, sizeof(size_t));
            output.write((char *)&val, sizeof(PRECISION));
        }
    }
}

SparseMatrix *binary_read_smat(std::istream &input)
{
    typedef Array2D_tuple<size_t, PRECISION> CCS_arr_type;

    size_t r, c, CCS_size, row_ind;
    PRECISION val;
    input.read((char *)&r, sizeof(size_t));
    input.read((char *)&c, sizeof(size_t));

    SparseMatrix *Ans = new SparseMatrix(r, c);
    Ans->initialize_CCS();
    for (size_t i = 0; i < c; i++)
    {
        input.read((char *)&CCS_size, sizeof(size_t));
        Array<size_t> *row_ind_i = new Array<size_t>(0, CCS_size);
        Array<PRECISION> *val_i = new Array<PRECISION>(0, CCS_size);
        for (size_t j = 0; j < CCS_size; j++)
        {
            input.read((char *)&row_ind, sizeof(size_t));
            input.read((char *)&val, sizeof(PRECISION));
            row_ind_i->push(row_ind);
            val_i->push(val);
        }
        Ans->push_CCS_col(row_ind_i, val_i);
    }
    return Ans;
}

void binary_print_lowertri(std::ostream &output, SpecMat::LowerTri<PRECISION> *ltmat)
{
    if (ltmat == nullptr)
    {
        unsigned dim = 0;
        output.write((char *)&dim, sizeof(unsigned));
        return;
    }

    unsigned dim = ltmat->dimension();
    output.write((char *)&dim, sizeof(unsigned));

    PRECISION *data = ltmat->get_vec();
    size_t ind = 0;
    PRECISION val;
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j <= i; j++)
        {
            val = data[ind++];
            output.write((char *)&val, sizeof(PRECISION));
        }
    }
}

SpecMat::LowerTri<PRECISION> *binary_read_lowertri(std::istream &input)
{
    size_t n;
    input.read((char *)&n, sizeof(size_t));
    if (n == 0)
        return nullptr;
    size_t size = (n * (n + 1)) / 2;
    PRECISION *val = new PRECISION[size];
    for (size_t i = 0; i < size; i++)
    {
        input.read((char *)&val[i], sizeof(PRECISION));
    }
    SpecMat::LowerTri<PRECISION> *Ans = new SpecMat::LowerTri(n, val);
    return Ans;
}

void print_lowertri(std::ostream &output, SpecMat::LowerTri<PRECISION> *ltmat)
{
    if (ltmat == nullptr)
        return;

    unsigned dim = ltmat->dimension();

    PRECISION *data = ltmat->get_vec();
    size_t ind = 0;
    PRECISION val;
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j <= i; j++)
        {
            output << data[ind++] << '\t';
        }
        output << '\n';
    }
}


void binary_print_taxon(std::ostream &output, TaxonList *taxa)
{
    output << taxa->size << taxa->bitstr_size;
    bool I2T, I2I;
    I2T = !taxa->Ind2Taxon.is_empty();
    I2I = !taxa->IndB2IndA.is_empty();

    int n = 0;

    output << I2T;
    if (I2T)
    {
        n = taxa->Ind2Taxon.get_size();
        output << n;
        for (size_t i = 0; i < n; i++)
            output << taxa->Ind2Taxon[i];
    }

    output << I2I;
    if (I2I)
    {
        n = taxa->IndB2IndA.get_size();
        output << n;
        for (size_t i = 0; i < n; i++)
            output << taxa->IndB2IndA[i];
    }
}

TaxonList *binary_read_taxon(std::istream &input)
{
    int b_size;
    input >> b_size; //droping the first size;
    input >> b_size;

    int t_size = 0, ind = 0;

    bool I2T, I2I;
    std::string temp;

    TaxonList *Ans = new TaxonList();
    Ans->set_bitstr_size(b_size);

    input >> I2T;
    if (I2T)
    {
        input >> t_size;
        for (int i = 0; i < t_size; i++)
        {
            input >> temp;
            Ans->push(temp, false);
        }
    }

    input >> I2I;
    if (I2I)
    {
        input >> t_size;
        for (int i = 0; i < t_size; i++)
        {
            input >> ind;
            Ans->IndB2IndA.push(ind);
        }
    }
    return Ans;
}
