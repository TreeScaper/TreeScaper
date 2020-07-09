
//##########################################################################
//# This software is part of the Treescaper i
//# -- Version 0.1
//# Copyright (C) 2010 Wen Huang
//#
//# This program is free software; you can redistribute it and/or
//# modify it under the terms of the GNU General Public License
//# as published by the Free Software Foundation; either version 2
//# of the License, or (at your option) any later version.
//#
//# This program is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//# http://www.gnu.org/copyleft/gpl.html
//##########################################################################

// TreeOPE.cpp
// Member function definitions for class TreeOPE.h
//              by whuang gzhou


//########################ZD comment########################################
//# This structure stores a tree by a set of linked nodes.
//# Most of the routines are implemented as recursive algorithms
//# in preorder traversal.
//########################ZD comment########################################

#ifndef TREEOPE_CPP
#define TREEOPE_CPP

#include "TreeOPE.h"

//########################ZD comment########################################
//# This is bit-wise addition and assignment macro-routines.
//########################ZD comment########################################

#define __HALF_MAX_SIGNED(type) ((type)1 << (sizeof(type)*8-2))
#define __MAX_SIGNED(type) (__HALF_MAX_SIGNED(type) - 1 + __HALF_MAX_SIGNED(type))
#define __MIN_SIGNED(type) (-1 - __MAX_SIGNED(type))
#define __MIN(type) ((type)-1 < 1?__MIN_SIGNED(type):(type)0)
#define __MAX(type) ((type)~__MIN(type))

//#ifndef _WIN32
    #define assign(dest,src) ({ \
    typeof(src) __x=(src); \
    typeof(dest) __y=__x; \
    (__x==__y && ((__x<1) == (__y<1)) ? (void)((dest)=__y),0:1); \
    })
    #define add_of(c,a,b) ({ \
    typeof(a) __a=a; \
    typeof(b) __b=b; \
    (__b)<1 ? \
    ((__MIN(typeof(c))-(__b)<=(__a)) ? assign(c,__a+__b):1) : \
    ((__MAX(typeof(c))-(__b)>=(__a)) ? assign(c,__a+__b):1); \
    })
//#else
//// This is not correct and need to be rewitten.
//    template<typename T1, typename T2, typename T3>
//    T3 assign(T1 dest, T2 src)
//    {
//        return dest;
//    }

//    template<typename T1, typename T2, typename T3>
//    T1 add_of(T1 c, T2 a, T3 b)
//    {
//        return a+b;
//    }

//    #define assign(dest,src) ({ \
//    typeid(src) __x=(src); \
//    typeid(dest) __y=__x; \
//    (__x==__y && ((__x<1) == (__y<1)) ? (void)((dest)=__y),0:1); \
//    })
//    #define add_of(c,a,b) ({ \
//    typeid(a) __a=a; \
//    typeid(b) __b=b; \
//    (__b)<1 ? \
//    ((__MIN(typeid(c))-(__b)<=(__a)) ? assign(c,__a+__b):1) : \
//    ((__MAX(typeid(c))-(__b)>=(__a)) ? assign(c,__a+__b):1); \
//    })
//#endif



/*
  load a  newick tree from a file
  Params: fname - path to file
          error - return for error code
  Returns: the loaded object on success, 0 on fail
  Error codes: -1 out of memory
               -2 parse error
               -3 can't open file
 */

 //########################ZD comment########################################
 //# Read newicktree from file "fname". The tree is stored
 //# as string, formatted by "()", same level a pair, ",",
 //# separator of same level of pairs and "\\" for quoted
 //# informations. Note that this routine only handles
 //# UNWEIGHTED tree.
 //########################ZD comment########################################

NEWICKTREE *TreeOPE::loadnewicktree(char *fname, int *error)
{
  NEWICKTREE *answer = NULL;
  FILE *fp = NULL;
  int err;

  fp = fopen(fname, "r");
  if(!fp)
  {
    err = -3;
    goto error_exit;
  }
  answer = floadnewicktree(fp, &err);
  if(!answer)
    goto error_exit;
  if(error)
    *error = 0;
  return answer;

 error_exit:
  if(error)
    *error = err;
  return 0;

}

//########################ZD comment########################################
//# This seems to be a duplicated loadnewicktree.
//# Both routines are not actually called in TreeScaper.
//########################ZD comment########################################

NEWICKTREE *TreeOPE::loadnewicktree2(FILE *fp, int *error)
{
  NEWICKTREE *answer = NULL;
//  FILE *fp;
  int err;

//  fp = fopen(fname, "r");
//  if(!fp)
//  {
//    err = -3;
//    goto error_exit;
//  }
  answer = floadnewicktree(fp, &err);
  if(!answer)
    goto error_exit;
  if(error)
    *error = 0;
  return answer;

 error_exit:
  if(error)
    *error = err;
  return 0;

}

/*
  load newick tree from an opened stream
  Params: fp - input stream
          error - return for error
  Returns: pointer to object, 0 on fail
  Error codes -1 out of memory
              -2 parse error
  Notes: format allows several trees to be stored in a file
 */

 //########################ZD comment########################################
 //# The implementation of the loadnewicktree. It process
 //# the string by calling subroutine loadnode in preorder
 //# recursively. This can be implemented by stack.
 //########################ZD comment########################################

NEWICKTREE *TreeOPE::floadnewicktree(FILE *fp, int *error)
{
  NEWICKTREE *answer = NULL;
  int err;
  int ch;

  answer = (NEWICKTREE *) malloc(sizeof(NEWICKTREE));
  if(!answer)
  {
    err = -1;
    goto error_exit;
  }
  skipspace(fp);
  ch = fgetc(fp);
  if(ch == '(' )
  {
    answer->root = loadnode(fp, &err);
    if(!answer->root || err != 0)
      goto error_exit;
  }
  skipspace(fp);
  ch = fgetc(fp);
  if(ch != ';')
  {
    err = -2;
    goto error_exit;
  }

  if(error)
    *error = 0;
  return answer;

 error_exit:
  killnewicktree(answer);
  if(error)
    *error = err;
  return 0;

}

/*
  newick tree destructor
 */

void TreeOPE::killnewicktree(NEWICKTREE *tree)
{

//    rewind(stdin);
   // setbuf(stdout, NULL);
  if(tree)
  {
    killnoder(tree->root);
    free(tree);
  }
}

/*
  node destructor (recursive)
 */

void TreeOPE::killnoder(NEWICKNODE *node)
{
  if(!node)
    return;

  for(int i = 0;i < node->Nchildren; i++)
  {
    killnoder(node->child[i]);
  }

  setbuf(stdout, NULL);

  free(node->label);
  if(node->bitstr != NULL)
  {
      delete node->bitstr;
      node->bitstr = NULL;
  }
  free(node->child);
  free(node);
}

void TreeOPE::killnewicktreeBitstr(NEWICKTREE *tree)
{
  if(tree)
  {
    killnoderBitstr(tree->root);
  }
}

/*
  node bitstr destructor (recursive)
 */

void TreeOPE::killnoderBitstr(NEWICKNODE *node)
{
  for(int i = 0;i < node->Nchildren; i++)
  {
      killnoderBitstr(node->child[i]);
  }

  if(node->bitstr != NULL)
  {
      delete node->bitstr;
      node->bitstr = NULL;
  }
}

/*
  turns a string into a label suitable for use in the tree
  Params: str - the string
  Returns: modified string, 0 on out of memory
  Notes: strings containing spaces have them replaced by underscores
         strings contianing illegal characters are quoted
         null pointer is returned as the empty string ""
 */
char *TreeOPE::makenewicklabel(const char *str)
{
  const char *cptr = NULL;
  char *vptr = NULL;
  int needsquote = 0;
  char *answer = NULL;

  if(!str)
    return mystrdup("");
  cptr = str;
  while(*cptr)
  {
    if(strchr("\'()[]:;,", *cptr))
      needsquote = 1;
    if(isspace(*cptr) && *cptr != ' ')
      needsquote = 1;
    cptr++;
  }
  if(needsquote)
  {
    answer = (char *) malloc(strlen(str) + 2 + mystrcount(str, '\''));
    vptr = answer;
    *vptr++ = '\'';
    while(*str)
    {
      *vptr = *str;
      if(*str == '\'')
      {
        vptr++;
        *vptr = '\'';
      }
      vptr++;
      str++;
    }
    *vptr++ = '\'';
    *vptr = 0;
  }
  else
  {
    answer = mystrdup(str);
    vptr = answer;
    while(*vptr)
    {
      if(*vptr == ' ' )
    *vptr = '_';
      vptr++;
    }
  }

  delete cptr;
  free(vptr);
  return answer;
}

/*
  load a node from the file
  Params: fp - the input stream
          error - return for error
  Returns: node loaded.
  Notes: recursive. Expects the opening parenthesis to have been eaten
 */

//########################ZD comment########################################
//# Recursive subroutine for floadnewicktree. "(" indicates
//# next level of pair, "," indicates the other node of the
//# current pair, ")" indicates returning to the previous
//# level. Note that when the string does not start with "(",
//# it means a leaf is reached and loadleaf should be called.
//########################ZD comment########################################

NEWICKNODE *TreeOPE::loadnode(FILE *fp, int *error)
{
  NEWICKNODE *answer = NULL;
  int err;
  NEWICKNODE *child = 0;
  int ch;

  answer = (NEWICKNODE *) malloc(sizeof(NEWICKNODE));
  //printf("sizeof(NEWICKNODE) = %d\n", sizeof(NEWICKNODE));
  if(!answer)
  {
    err = -1;
    goto error_exit;
  }

  answer->Nchildren = 0;
  answer->label = 0;
  answer->child = 0;
  answer->hv1 = 0;
  answer->hv2 = 0;
  answer->bitstr = NULL;

  skipspace(fp);
  do
  {
    ch = fgetc(fp);
    if(ch == '(')
    {
      child = loadnode(fp, &err);
      if(!child)
        goto error_exit;

      if(addchild(answer, child ) == -1)
      {
        err = -1;
        goto error_exit;
      }
      child = 0;
    }
    else
    {
      ungetc(ch, fp);
      child = loadleaf(fp, &err);
      if(!child)
        goto error_exit;

      if(addchild(answer, child) == -1)
      {
        err = -1;
        goto error_exit;
      }
      child = 0;
    }
    skipspace(fp);
    ch = fgetc(fp);
  } while(ch == ',');

  if(ch == ')')
  {
    answer->label = readlabelandweight(fp, &answer->weight, &err);
    if(err)
      goto error_exit;
  }
  else
  {
    err = -2;
    goto error_exit;
  }

  if(error)
    *error = 0;
  return answer;

 error_exit:
  if(child)
    killnoder(child);
  killnoder(answer);
  if(error)
    *error = err;
  return 0;
}

/*
  add a child to a node
  Params: parent - parent node
          child - child to add
  Returns: 0 on success, -1 on fail
 */
int TreeOPE::addchild(NEWICKNODE *parent, NEWICKNODE *child)
{
  NEWICKNODE **temp = NULL;

  temp = (NEWICKNODE **) realloc(parent->child, (parent->Nchildren + 1) * sizeof(NEWICKNODE *));
//  std::cout << "addchild:" << (long) temp << std::endl;//-- WHtest
//  if(child->label != '\0')
//      std::cout << "cname:" << child->label << std::endl;//--- WHtest
  if(!temp)
    return -1;
  parent->child = temp;
  parent->child[parent->Nchildren] = child;
  parent->Nchildren++;

  return 0;
}


/*
  load a leaf node
  Params: fp - the input stream
          error - return for error
  Returns: node object

 */
NEWICKNODE *TreeOPE::loadleaf(FILE *fp, int *error)
{
  NEWICKNODE *answer = NULL;
  int err;

  answer = (NEWICKNODE *) malloc(sizeof(NEWICKNODE));
  if(!answer)
  {
    err = -1;
    goto error_exit;
  }

  answer->Nchildren = 0;
  answer->child = 0;
  answer->label = readlabelandweight(fp, &answer->weight, &err);
  answer->bitstr = NULL;
  if(err)
    goto error_exit;
  return answer;

  error_exit:
  if(error)
    *error = err;
  if(answer)
  {
    free(answer->label);
    free(answer);
  }
  return 0;
}

/*
  read label, colon and associated weight
  Params: fp - the input stream
          weight - return for weight
          error - return for error
  Returns: the label
  Notes: a null label is not an error
 */
char *TreeOPE::readlabelandweight(FILE *fp, double *weight, int *error)
{
  char *answer = NULL;
  int ch;

  *weight = 0;
  answer = readlabel(fp);
  if(!answer)
  {
    *error = -1;
    return 0;
  }
  skipspace(fp);
  ch = fgetc(fp);
  if(ch == ':')
  {
    if( fscanf(fp, "%lf", weight) != 1)
    {
      *error = -2;
      free(answer);
      return 0;
    }
  }
  else
    ungetc(ch, fp);

  if(*answer == 0)
  {
    free(answer);
    answer = 0;
  }
  *error = 0;
  return answer;
}

/*
  readlabel - read a label from stream
  Params: fp - the input stream
  Returns: allocarted label
  Notes: null label is the null string ""
 */
char *TreeOPE::readlabel(FILE *fp)
{
  char *answer = NULL;
  int ch;

  skipspace(fp);
  ch = fgetc(fp);
  ungetc(ch, fp);
  if(ch == '\'')
    answer = readquotedlabel(fp);
  else
    answer = readplainlabel(fp);

  return answer;
}

/*
  read a quoted label from stream
  Params: fp - the input stream
  Returns: the label. null string is ""
 */
char *TreeOPE::readquotedlabel(FILE *fp)
{
  char *answer = NULL;
  char *temp = NULL;
  int capacity = 32 * sizeof(char);
  int len = 0;
  int ch;

  answer = (char *) malloc(capacity);
  if(!answer)
    return 0;
  ch = fgetc(fp);
  while( (ch = fgetc(fp)) != EOF)
  {
    if(ch == '\'')
    {
      ch = fgetc(fp);
      if(ch == '\'')
        answer[len++] = (char) ch;
      else
      {
        ungetc(ch, fp);
        answer[len] = 0;
        return answer;
      }
    }
    else
      answer[len++] = (char) ch;
    if(len == capacity - 1)
    {
      temp = (char *) realloc(answer, capacity * 2);
      if(!temp)
      {
        free(answer);
        return 0;
      }
      answer = temp;
      capacity = capacity * 2;
    }
  }
  answer[len] = 0;
  return answer;
}

/*
  read aplian (unquoted) label
  Params: fp - the open file
  Returns: the label, null string is ""
 */
char *TreeOPE::readplainlabel(FILE *fp)
{
  char *answer = NULL;
  char *temp = NULL;
  int len = 0;
  int capacity = 32 * sizeof(char);
  int ch;

  answer = (char *) malloc(capacity);
  if(!answer)
    return 0;

  while( (ch = fgetc(fp)) != EOF)
  {
    if(isspace(ch) || strchr(" ()[]\':;,", ch))
    {
      ungetc(ch, fp);
      answer[len] = 0;
      return answer;
    }
    if(ch == '_')
      ch = '_'; // ssj
    answer[len++] = (char) ch;
    if(len == capacity - 1)
    {
      temp = (char *) realloc(answer, capacity * 2);
      if(!temp)
      {
        free(answer);
        return 0;
      }
      answer = temp;
      capacity = capacity * 2;
    }
  }

  answer[len] = 0;
  return answer;
}

/*
  consume space in input stream.
  Params: fp - the input stream
 */
void TreeOPE::skipspace(FILE *fp)
{
  int ch;

  while( (ch = fgetc(fp)) != EOF )
    if(!isspace(ch))
    {
      ungetc(ch, fp);
      break;
    }
}

/*
  strduo implemetation
 */
char *TreeOPE::mystrdup(const char *str)
{
  char *answer = NULL;

  answer = (char *) malloc(strlen(str) + 1);
  if(answer)
    strcpy(answer, str);

  return answer;
}

/*
  count the number of instances of ch in the string
  Params: str - serach string
          ch - character.
  Returns: number of times ch appears in str.
 */
int TreeOPE::mystrcount(const char *str, int ch)
{
  int answer = 0;

  while(*str)
    if(*str++ == ch)
      answer++;

  return answer;
}

/*
  print out a node (test function)
  Parmas: node - node to print
          indent - indent level
  Notes: recursive
*/
void TreeOPE::printnewicknode(NEWICKNODE *node, int indent)
{
  int i;
  for(i=0;i<indent;i++)
      cout << " ";
  cout << (long) node << ":";
  if(node->label == NULL)
      cout << "NULL";
  else
      cout << node->label;
  cout << " : " << node->weight << " : " << node->Nchildren << endl;

  for(i=0;i<node->Nchildren;i++)
    printnewicknode(node->child[i], indent + 5);
}

void TreeOPE::printnewicktree(const NEWICKTREE *tree)
{
  printnewicknode(tree->root, 0);
}

int TreeOPE::newickmain(int argc, char **argv)
{
  NEWICKTREE *tree = NULL;
  int err;
  char *label = NULL;

  if(argc != 2)
    fprintf(stderr, "Newick tree loader\n");
  else
  {
    tree = loadnewicktree(argv[1], &err);
    if(!tree)
    {
      switch(err)
      {
      case -1: printf("Out of memory\n"); break;
      case -2: printf("parse error\n"); break;
      case -3: printf("Can't load file\n"); break;
      default:
        printf("Error %d\n", err);
      }
    }
    else
    {
      printf("Loaded\n");
      printnewicktree(tree);
    }
    killnewicktree(tree);
  }
  label = makenewicklabel("My name is \'Fred()\'");
  printf("***%s***\n", label);
  free(label);

  return 0;
}

void TreeOPE::GetTaxaLabels(NEWICKNODE *node, LabelMap &lm)
{
    if (node->Nchildren == 0)
    {
        std::string temp(node->label);
        lm.push(temp);
    }
    else
        for (int i = 0; i<node->Nchildren; i++)
            GetTaxaLabels(node->child[i], lm);
};

//########################ZD comment########################################
//# Assigned hash values of computed from its sub-trees,
//# in recursive way. For leaves, hv1 and hv2 are read
//# from the hash table directly. For internal nodes,
//# the hash values are the sum of its children's hash
//# values. The implementation here is unnecessarily
//# complicated. Also note that the sum of hash values
//# is done by the bit-wise marco addition defined earlier.
//# The particular hash value hv2 is associated to
//# this string representation of the bipartition and
//# stored in hash2bitstr.
//########################ZD comment########################################

void TreeOPE::dfs_compute_hash(
                            NEWICKNODE* startNode,
                            LabelMap &lm,
                            HashRFMap &vec_hashrf,
                            unsigned treeIdx,
                            unsigned &numBitstr,
                            unsigned long long m1,
                            unsigned long long m2,
                            bool WEIGHTED,
                            unsigned int NUM_Taxa,
                            map<unsigned long long, Array<char> *> &hash2bitstr,
                            int numofbipartions)
{
    // If the node is leaf node, just set the place of the taxon name in the bit string to '1'
    // push the bit string into stack

    /*cout << "output test========" << endl;//---
    for(int i = 0; i < 16; i++)
    {
        Array<char> *arr = new Array<char> (16);
        arr->SetBitArray(i+1);
        arr->printbits(16);
        cout << endl;//---
        delete arr;
    }*/

    int btlength = (int) ((NUM_Taxa + 7) / 8);
    if (startNode->Nchildren == 0)
    {
        // leaf node
//        string temp(startNode->label);
//        unsigned int idx = lm[temp];
        unsigned int idx = atoi(startNode->label) - 1;

        // Implicit BPS
        // Set the hash value for each leaf node
        startNode->hv1 = vec_hashrf._HF.getA1(idx);
        startNode->hv2 = vec_hashrf._HF.getA2(idx);

//        std::cout << "l:" << startNode->bitstr.get_length() << ", idx:" << idx << std::endl;//--
        startNode->bitstr = new Array<char> (btlength);

//        for(int i = 0; i < btlength; i++)
//        {
//            std::cout << "s i:" << i << std::endl;//--
//            startNode->bitstr[i] = NULL;
//        }
//        cout << "ith:" << idx + 1 << endl;//----
        if(! startNode->bitstr->SetBitArray(idx + 1))
            std::cout << "Warning: Errors in setting bit string" << std::endl;
//        startNode->bitstr->printbits(12);//--
//        cout << endl;//---
    }
    else
    {
        for (int i = 0; i < startNode->Nchildren; ++i)
        {
            dfs_compute_hash(startNode->child[i], lm, vec_hashrf,
                             treeIdx, numBitstr, m1, m2,WEIGHTED,NUM_Taxa, hash2bitstr, numofbipartions);
        }
        // For weighted RF
        float dist = 0.0;
        if (WEIGHTED)
            dist = startNode->weight;
        else
            dist = 1;

        ++ numBitstr;

        // Implicit BPS
        // After an internal node is found, compute hv1 and hv2
        unsigned long long temphv1 = 0;
        unsigned long long temphv2 = 0;
        startNode->bitstr = new Array<char> (btlength);

        for (int i = 0; i < startNode->Nchildren; ++i)
        {
//            cout << "i:" << i << ",child:" << endl;///---
//            startNode->child[i]->bitstr->printbits(12);//---
            Array<char> bitstrleft = *(startNode->bitstr);
            Array<char> bitstrright = *(startNode->child[i]->bitstr);
            startNode->bitstr->ORbitOPE(bitstrleft, bitstrright);
//            cout << endl;//--
//            cout << "results:";
//            startNode->bitstr->printbits(12);//----
//            cout << endl;//----
            unsigned long long t1 = temphv1;
            unsigned long long t2 = temphv2;
            unsigned long long h1 = startNode->child[i]->hv1;
            unsigned long long h2 = startNode->child[i]->hv2;

            if (add_of(temphv1, t1, h1))
            {
                cout << "ERROR : ullong add overflow !\n\n";
                cout << "t1 = " << t1 << "h1 = " << h1 << "t1 + h1 = " << t1+h1 << endl;
                exit(0);
            }
            if (add_of(temphv2, t2, h2))
            {
                cout << "ERROR : ullong add overflow !\n";
                cout << "t2 = " << t2 << "h2 = " << h2 << "t2 + h2 = " << t2+h2 << endl;
                exit(0);
            }
        }
//        startNode->bitstr->printbits(12);//---
//        cout << endl;//---

        // Check overflow
        unsigned long long temp1 = temphv1 % m1;
        unsigned long long temp2 = temphv2 % m2;
        startNode->hv1 = temp1;
        startNode->hv2 = temp2;

        Array<char> *btpt = NULL;
        btpt = new Array<char> (*startNode->bitstr);
        if(hash2bitstr[startNode->hv2] != NULL)
        {
            delete hash2bitstr[startNode->hv2];
            hash2bitstr[startNode->hv2] = btpt;
        }
        else
        {
            hash2bitstr[startNode->hv2] = btpt;
        }
        // Store bitstring in hash table
        if (numBitstr <= numofbipartions)
        {
            vec_hashrf.hashing_bs_without_type2_nbits(treeIdx, NUM_Taxa, startNode->hv1, startNode->hv2, dist, WEIGHTED);
        }
    }
};

//########################ZD comment########################################
//# Create the arrays that store hash values,
//# TreeIdx and weight in preorder. Note that
//# the "TreeIdx" is an identical array.
//# Each tree will generate one set of such
//# arrays and arrays from different trees
//# are pasted together and sorted by the hash
//# values. By comparing hash values, identical
//# bipartitions are found.
//########################ZD comment########################################

void TreeOPE::bipart(NEWICKNODE * const startnode, unsigned int &treeIdx,
                         unsigned long long *matrix_hv,
                         unsigned int *matrix_treeIdx,
                         double *matrix_weight, int &idx, int depth, bool isrooted)
{
    depth++;
    if (isrooted)
    {
        if(depth > 1 && startnode->Nchildren > 0)
        {
            matrix_hv[idx] = startnode->hv2;
            matrix_treeIdx[idx] = treeIdx+1;
            matrix_weight[idx] = startnode->weight;
            idx++;
        }
    }
    else
    {   if(depth > 2 && startnode->Nchildren > 0)
        {
            matrix_hv[idx] = startnode->hv2;
            matrix_treeIdx[idx] = treeIdx+1;
            matrix_weight[idx] = startnode->weight;
            idx++;
        }
    }

    for (int i = 0; i < startnode->Nchildren; i++)
        if(startnode->child[i]->Nchildren > 0)
            bipart(startnode->child[i], treeIdx, matrix_hv, matrix_treeIdx, matrix_weight, idx, depth, isrooted);
};

//########################ZD comment########################################
//# This is a duplicated code for reading weighted binary
//# tree, which is actually called by TreeScaper.
//########################ZD comment########################################

NEWICKTREE *TreeOPE::parsetree(char *str, int *error, NEWICKTREE *testtree)
{
    NEWICKTREE *answer = NULL;
    int err;
    char *ch = NULL;
    int i = 0, idx = 0;

    answer = (NEWICKTREE *) malloc(sizeof(NEWICKTREE));

    // skip space
    char *str_copy = (char *) malloc(sizeof(char) * (strlen(str) + 2));
    for(i = 0; i < strlen(str) + 1; i++)
    {
        if((*(str + i)) != ' ')
        {
            *(str_copy + idx) = *(str + i);
            idx++;
        }
        if(*(str + i) == '\0')
            break;
    }
    *(str_copy + idx) = '\0';

    if (!answer)
    {
        err = -1;
        goto error_exit;
    }

    ch = str_copy;
    if (*ch == '(')
    {
        ch++;
        answer->root = parsenode(&ch, &err);
        if (!answer->root || err != 0)
            goto error_exit;
    }

    if(testtree != NULL)
    {
        cout << "in2-----" << endl;
        TreeOPE::printnewicktree(testtree);
    }

    if (*ch != ';')
    {
        err = -2; // parse error
        goto error_exit;
    }

    if (error)
        *error = 0;
    free(str_copy);
//    std::cout << "root->child:" << (long) answer->root->child << std::endl;//-- WHtest
    return answer;

error_exit:
    cout << "Wrong in parsetree()!\n\n";//----
    free(str_copy);
    killnewicktree(answer);

    if (error)
        *error = err;
    return NULL;
}

NEWICKNODE *TreeOPE::parsenode(char **str, int *error)
{
    NEWICKNODE *answer = NULL;
    int err;
    NEWICKNODE *child = NULL;

    answer = (NEWICKNODE *) malloc(sizeof(NEWICKNODE));
    if (!answer)
    {
        err = -1;
        goto error_exit;
    }

    answer->Nchildren = 0;
    answer->label = NULL;
    answer->child = NULL;
    answer->hv1 = 0;
    answer->hv2 = 0;
    answer->bitstr = NULL;

    if (**str == '(')
    {
        (*str)++;
        child = parsenode(str, &err);
        if (!child)
            goto error_exit;

        if (addchild(answer, child) == -1)
        {
            err = -1;
            goto error_exit;
        }
        child = NULL;
    }
    else
    {
        child  = parseleaf(str, &err);
        if (!child)
            goto error_exit;

        if (addchild(answer, child) == -1)
        {
            err = -1;
            goto error_exit;
        }
        child = NULL;
    }

    while(**str ==  ',')
    {
        (*str)++;
        if (**str == '(')
        {
            (*str)++;
            child = parsenode(str, &err);

            if (!child)
                goto error_exit;

            if (addchild(answer, child) == -1)
            {
                err = -1;
                goto error_exit;
            }
            child = NULL;
        }
        else
        {
            child  = parseleaf(str, &err);
            if (!child)
                goto error_exit;

            if (addchild(answer, child) == -1)
            {
                err = -1;
                goto error_exit;
            }
            child = NULL;
        }
    }

    if (**str == ')')
    {
        (*str)++;
        answer->label = parselabelandweight(str, &answer->weight, &err);
        if(err)
            goto error_exit;
    }
    else
    {
        err = -2;
        goto error_exit;
    }

    if (error)
        *error = 0;
    return answer;

error_exit:
    cout << "Wrong in parsenode()!\n\n";//----
    if (child)
        killnoder(child);
    killnoder(answer);

    if (error)
        *error = err;
    return NULL;
}

NEWICKNODE *TreeOPE::parseleaf(char **str, int *error)
{
    NEWICKNODE *answer = NULL;
    int err;

    answer = (NEWICKNODE *) malloc(sizeof(NEWICKNODE));
    if (!answer)
    {
        err = -1;
        goto error_exit;
    }
    answer->Nchildren = 0;
    answer->child = NULL;

    answer->label = parselabelandweight(str, &answer->weight, &err);
    answer->bitstr = NULL;

    if (err)
        goto error_exit;
    return answer;

error_exit:
    cout << "Wrong in parseleaves()!\n\n";//----
    if (error)
        *error = err;
    if (answer)
    {
        free(answer->label);
        free(answer);
    }
    return NULL;
}

char *TreeOPE::parselabelandweight(char **str, double *weight, int *error)
{
    char *answer = NULL;
    int a;
    *weight = 0;
//    cout << "str:" << *str << endl;//----
    answer = parselabel(str);

//    cout << "str:" << *str << endl;//----
//    cout << "label:" << answer << endl;//----
    if (!answer)
    {
        cout << "Wrong in parselabelandweight()!\n\n";//----
        *error = -1;
        return 0;
    }

    if (**str == ':')
    {
        (*str)++;
        a = sscanf(*str, "%lf", weight);
        while ((**str) == '.' || ((**str) >= '0' && (**str) <= '9') || (**str) == 'e' || (**str) == '-' || (**str) == 'E' || (**str) == '+')
            (*str)++;
        if (a != 1)
        {
            cout << "Wrong in parselabelandweight()!\n\n";//----
            *error = -2;
            free(answer);
            return 0;
        }
    }

    if (*answer == NULL)
    {
        free(answer);
        answer = NULL;
    }
    *error = 0;
    return answer;
}

char *TreeOPE::parselabel(char **str)
{
    char *answer = NULL;
    if (**str == '\'')
    {
        answer = parsequotedlabel(str);
    }
    else
    {
        answer = parseplainlabel(str);
    }

    return answer;
}

char *TreeOPE::parsequotedlabel(char **str)
{
    char *answer = NULL;
    char *temp = NULL;
    int capacity = 32 * sizeof(char);
    int len = 0;

    answer = (char *) malloc(capacity);
    if (!answer)
    {
        cout << "Wrong in parsequotedlabel()!\n\n";//----
        return 0;
    }

    (*str)++;
    while (**str != '\0')
    {
        if (**str == '\'')
        {
            (*str)++;
            if (**str == '\'')
            {
                answer[len++] = **str;
            }
            else
            {
                answer[len] = 0;
                return answer;
            }
        }
        else
            answer[len++] = **str;
        if (len == capacity - 1)
        {
            temp = (char *) realloc(answer, capacity * 2);
            if (!temp)
            {
                free(answer);
                cout << "Wrong in parsequotedlabel()!\n\n";//----
                return 0;
            }
            answer = temp;
            capacity = capacity * 2;
        }
    }
    answer[len] = 0;
    return answer;
}

char *TreeOPE::parseplainlabel(char **str)
{
    char *answer = NULL;
    char *temp = NULL;
    int len = 0;
    int capacity = 32 * sizeof(char);

    answer = (char *) malloc(capacity);
    if (!answer)
    {
        cout << "Wrong in parseplainlabel()!\n\n";//----
        return 0;
    }

    while (**str != '\0')
    {
        if (isspace((unsigned char)**str) || strchr(" ()[]\':;,", (**str)))
//        if((**str) == ' ' || (**str) == '(' || (**str) == ')' || (**str) == '[' || (**str) == ']'
//                || (**str) == '\'' || (**str) == ':' || (**str) == ';' || (**str) == ',')
        {
            answer[len] = 0;
            return answer;
        }

        if (**str == '_')
            **str = '_';
        answer[len++] = **str;
        if (len == capacity - 1)
        {
            temp = (char *) realloc(answer, capacity * 2);
            if (!temp)
            {
                free(answer);
                return 0;
            }
            answer = temp;
            capacity = capacity * 2;
        }
        (*str)++;
    }
    answer[len] = 0;
    return answer;
}

/*
 print the Newick tree in a human-readable format (test function)
 Params: tree - the tree

void TreeOPE::printnewicknode(NEWICKNODE *node, int indent)
{
    int i;
    if (node->Nchildren == 0)
        printf("%s:%f",node->label, node->weight);
    else
    {
        printf("(");
        for (i = 0; i < node->Nchildren; i++)
        {
            printnewicknode(node->child[i], indent);
            if (i != node->Nchildren-1)
                printf(",");
            else
                printf(")");
        }
    }
}

void TreeOPE::printnewicktree(NEWICKTREE *tree)
{
    printnewicknode(tree->root, 0);
}
*/
/* print Newick tree in Nexus format*/
void TreeOPE::printTree_nex(NEWICKNODE *node, int length, char *gzbuff, int depth, int N)
{
    depth++;
    int bufflen, i;
    if (node->Nchildren == 0)
    {
        bufflen = strlen(gzbuff);
        if (bufflen > N - 100)
        {
            gzbuff = (char*)realloc(gzbuff, 2 * N);
        }
        if (node->weight == 0.0)
            sprintf(gzbuff + bufflen, "%s:%0.6f", node->label, 1.0);
        else
            sprintf(gzbuff + bufflen, "%s:%0.6f", node->label, node->weight);
    }
    else
    {
        bufflen = strlen(gzbuff);
        sprintf(gzbuff + bufflen, "(");
        for (i = 0; i < node->Nchildren; i++)
        {
            printTree_nex(node->child[i], length, gzbuff, depth, N);
            if (i != node->Nchildren-1)
            {
                bufflen = strlen(gzbuff);
                if (bufflen > N - 100)
                    gzbuff = (char*) realloc(gzbuff, 2 * N);
                sprintf(gzbuff + bufflen, ",");
            }
        }

        if (node->label != NULL)
        {
             bufflen = strlen(gzbuff);
             if (bufflen > N - 100)
                 gzbuff = (char*)realloc(gzbuff, 2 * N);
             if (depth != 1)
             {
                 if (node->weight == 0.0)
                     sprintf(gzbuff + bufflen, "%s:%0.6f", node->label, 1.0);
                 else
                     sprintf(gzbuff + bufflen, "%s:%0.6f", node->label, node->weight);
             }
             else
                 sprintf(gzbuff + bufflen, "%s", node->label);
        }
        else
        {
            bufflen = strlen(gzbuff);
            if (bufflen > N - 100)
                gzbuff = (char*)realloc(gzbuff, 2 * N);
            if (depth != 1)
            {
                if (node->weight == 0.0)
                    sprintf(gzbuff + bufflen, "):%0.7f", 1.0);
                else
                    sprintf(gzbuff + bufflen, "):%0.7f", node->weight);
            }
            else
                sprintf(gzbuff + bufflen, ")");

        }
    }
}

/* print Nexus tree in Newick format*/
void TreeOPE::printTree_new(NEWICKNODE *node, char **str, char *gzbuff, int depth, int N)
{
    depth++;
    int bufflen, i;
    if (node->Nchildren == 0)
    {
        bufflen = strlen(gzbuff);
        if (bufflen > N - 100)
            gzbuff = (char*)realloc(gzbuff, 2*N);
        if (node->weight == 0.0)
            sprintf(gzbuff + bufflen, "%s:%0.6f", str[atoi(node->label)-1], 1.0);
        else
            sprintf(gzbuff + bufflen, "%s:%0.6f", str[atoi(node->label)-1], node->weight);
    }
    else
    {
        bufflen = strlen(gzbuff);
        sprintf(gzbuff + bufflen, "(");
        for (i = 0; i < node->Nchildren; i++)
        {
            printTree_new(node->child[i], str, gzbuff, depth, N);
            if (i != node->Nchildren-1)
            {
                bufflen = strlen(gzbuff);
                if (bufflen > N - 100)
                    gzbuff = (char*) realloc(gzbuff, 2 * N);
                sprintf(gzbuff + bufflen, ",");
            }
        }
        if (node->label != NULL)
        {
            bufflen = strlen(gzbuff);
            if (bufflen > N - 100)
                gzbuff = (char*)realloc(gzbuff, 2*N);
            if (depth != 1)
            {
                if (node->weight == 0.0)
                    sprintf(gzbuff + bufflen, "%s:%0.6f", str[atoi(node->label)-1], 1.0);
                else
                    sprintf(gzbuff + bufflen, "%s:%0.6f", str[atoi(node->label)-1], node->weight);
            }
            else
                sprintf(gzbuff + bufflen, "%s", str[atoi(node->label) - 1]);
        }
        else
        {
            bufflen = strlen(gzbuff);
            if (bufflen > N-100)
            {
                gzbuff = (char*)realloc(gzbuff, 2*N);
            }
            if(depth != 1)
            {
                if (node->weight == 0.0)
                    sprintf(gzbuff + bufflen, "):%0.6f", 1.0);
                else
                    sprintf(gzbuff + bufflen, "):%0.6f", node->weight);
            }
            else
                sprintf(gzbuff + bufflen, ")");
        }
    }
}


NEWICKNODE *TreeOPE::findleaf(std::string leafname, NEWICKNODE *currentnode, NEWICKNODE *parent, int *icpt)
{
    currentnode->parent = parent;
    NEWICKNODE *resultnode = NULL;
    NEWICKNODE *temp = NULL;

    for(int i = 0; i < currentnode->Nchildren; i++)
    {
        temp = findleaf(leafname, currentnode->child[i], currentnode, icpt);
        if(temp != NULL)
        {
            resultnode = temp;
            if(parent == NULL)
                (*icpt) = i;
        }
    }

    if(currentnode->Nchildren == 0)
        if(0 == strcmp(leafname.c_str(), currentnode->label))
            resultnode = currentnode;

    return resultnode;
};

//########################ZD comment########################################
//# Lift the unrooted tree to rooted tree.  When a tree
//# is lifted, the node that become parent node does not
//# free from the memory.
//########################ZD comment########################################


void TreeOPE::normailzedTree(NEWICKNODE *lrpt, NEWICKTREE *newickTree, int indexchild)
{
    if(newickTree->root->Nchildren == 2)
    {
        newickTree->root->child[1 - indexchild]->parent = newickTree->root->child[indexchild];
        newickTree->root->child[1 - indexchild]->weight += newickTree->root->child[indexchild]->weight;
        newickTree->root->child[indexchild]->parent = newickTree->root->child[1 - indexchild];
        free(newickTree->root->child);
        newickTree->root->child = NULL;
        free(newickTree->root->label);
        free(newickTree->root);
    }

    NEWICKNODE *root = (NEWICKNODE *) malloc(sizeof(NEWICKNODE));
    root->child = (NEWICKNODE **) malloc(sizeof(NEWICKNODE *) * 2);
    root->child[0] = lrpt;
    root->child[1] = lrpt->parent;
    root->Nchildren = 2;
    root->parent = NULL;
    root->label = NULL;
    root->hv1 = 0;
    root->hv2 = 0;
    root->weight = 0.0;
    root->bitstr = NULL;
    newickTree->root = root;

    for(int i = 0; i < lrpt->parent->Nchildren; i++)
    {
        if(lrpt->parent->child[i] == lrpt)
        {
            lrpt->parent->child[i] = root;
            break;
        }
    }
    lrpt->parent = root;

    normalizedNode(newickTree->root, NULL, 0);
}

//########################ZD comment########################################
//# Recursive step of normailzedTree.
//########################ZD comment########################################

void TreeOPE::normalizedNode(NEWICKNODE *currentnode, NEWICKNODE *parent, double currentweight)
{
    NEWICKNODE *temp = NULL;
    double childweight = 0;

    if(currentnode->parent != parent)
    {
        if(currentnode->parent == NULL)
        {
            currentnode->parent = parent;
            for(int i = 0; i < currentnode->Nchildren; i++)
            {
                if(currentnode->child[i] == parent)
                {
                    currentnode->child[i] = currentnode->child[currentnode->Nchildren - 1];
                    break;
                }
            }
            currentnode->Nchildren = currentnode->Nchildren - 1;
            currentnode->weight = currentweight;
        } else
        {
            temp = currentnode->parent;
            currentnode->parent = parent;
            for(int i = 0; i < currentnode->Nchildren; i++)
            {
                if(currentnode->child[i] == parent)
                {
                    currentnode->child[i] = temp;
                    break;
                }
            }
            childweight = currentnode->weight;
            currentnode->weight = currentweight;
        }
    }
    for(int i = 0; i < currentnode->Nchildren; i++)
    {
        normalizedNode(currentnode->child[i], currentnode, childweight);
    }
}

void TreeOPE::Label_strint(NEWICKNODE *node, LabelMap &lm)
{
    if(node->Nchildren == 0)
    {
        std::ostringstream ostr;
        ostr << (int) lm[node->label] + 1;
        std:: string intstr = ostr.str();
        free(node->label);
        node->label = new char [intstr.size() + 1];
        strcpy(node->label, intstr.c_str());
    }

    for(int i = 0; i < node->Nchildren; i++)
    {
        Label_strint(node->child[i], lm);
    }
}

//########################ZD comment########################################
//# Convert Newicktree to Ptree, which is a index-
//# based tree structure, which does not stored the
//# the hash value and weight, i.e., the bipartition
//# and weight information are lost. This is used to
//# compute match distance. Recall that Ptree has
//# 3 arrays, indices of left child, indices of right
//# child and indices of parent and edges matrix,
//# which is not computed here.
//########################ZD comment########################################

bool TreeOPE::newick2lcbb(const NEWICKTREE *nwtree, int num_leaves, struct Ptree *tree)
{
    tree->leaf_number = num_leaves;
    for(int i = 0; i < 2 * num_leaves; i++)
    {
        tree->lchild[i] = -1;
        tree->rchild[i] = -1;
        for(int j = 0; j < 2 * num_leaves; j++)
            tree->edge[i][j] = - 1;
    }

    bool currentnode;
    int currentidx = -1;
    int vidx = num_leaves + 1;

    if (nwtree->root->Nchildren == 3)
    {
        NEWICKNODE *pseudo_root = new NEWICKNODE;
        NEWICKNODE *temp = new NEWICKNODE;
        pseudo_root->Nchildren = 2;
        pseudo_root->child = new NEWICKNODE *[2];
        pseudo_root->child[1] = temp;
        pseudo_root->child[0] = nwtree->root->child[0];
        temp->Nchildren = 2;
        temp->child = new NEWICKNODE *[2];
        temp->child[0] = nwtree->root->child[1];
        temp->child[1] = nwtree->root->child[2];

        if( !newick2ptree(pseudo_root, tree, currentnode, vidx, currentidx))
            return false;
        tree->parent[num_leaves] = -1;
        tree->lchild[num_leaves] = tree->lchild[2 * num_leaves - 1];
        tree->rchild[num_leaves] = tree->rchild[2 * num_leaves - 1];

        delete [] pseudo_root->child;
        delete pseudo_root;
        delete [] temp->child;
        delete temp;

        return true;
    }

    if(!newick2ptree(nwtree->root, tree, currentnode, vidx, currentidx))
        return false;

    tree->parent[num_leaves] = -1;
    tree->lchild[num_leaves] = tree->lchild[2 * num_leaves - 1];
    tree->rchild[num_leaves] = tree->rchild[2 * num_leaves - 1];
    return true;
}

//########################ZD comment########################################
//# Implementation of newick2lcbb. Leaves are
//# labeled with the additional information, taxa#
//# I assumed, and then tree is reconstructed in
//# array structure that stored left, right child
//# and parent label. Note that the edge matrix
//# is not updated here.
//########################ZD comment########################################

bool TreeOPE::newick2ptree(NEWICKNODE *node, struct Ptree *tree, bool &currentnode, int &vidx, int &curidx)
{
    if(node->Nchildren != 2 && node->Nchildren != 0)
    {
        return false;
    }
    bool child;
    currentnode = true;
    int cidx;

    int *childrenidx = new int [node->Nchildren];

    for(int i = 0; i < node->Nchildren; i++)
    {
        newick2ptree(node->child[i], tree, child, vidx, childrenidx[i]);
        currentnode = currentnode && child;
    }

    if(node->Nchildren == 0)
    {
        cidx = atoi(node->label) - 1;
        tree->lchild[cidx] = -1;
        tree->rchild[cidx] = -1;
        curidx = cidx;
        return true;
    }

    if(currentnode)
    {
        tree->lchild[vidx] = childrenidx[0];
        tree->rchild[vidx] = childrenidx[1];
        curidx = vidx;
        if(childrenidx[0] < 0 || childrenidx[1] < 0 || childrenidx[0] > vidx || childrenidx[1] > vidx)
        {
            return false;
        }

        tree->parent[childrenidx[0]] = curidx;
        tree->parent[childrenidx[1]] = curidx;

        vidx++;
        return true;
    }
    delete [] childrenidx;
    return true;
}

int TreeOPE::sumofdegree(NEWICKNODE *node, bool isrooted, int depth)
{
    depth++;
    int result = 0;
    if(depth == 0 && isrooted)
        result += node->Nchildren;
    if(depth > 0)
        result += node->Nchildren + 1;
    for(int i = 0; i < node->Nchildren; i++)
        result += sumofdegree(node->child[i], isrooted, depth);

    return result;
}

//########################ZD comment########################################
//# Redundant function
//########################ZD comment########################################

void TreeOPE::obtainbipartcount(NEWICKTREE *nwtree, bool isrooted, map<unsigned long long, unsigned long long> &bipcount)
{

    bipartcount(nwtree->root, isrooted, bipcount);
}

//########################ZD comment########################################
//# Use hash values to count bipartition. Recall that
//# hv2 represents a bipartition, accumulate the occurance
//# each hv2 and stored it in bipcount.
//########################ZD comment########################################

void TreeOPE::bipartcount(NEWICKNODE *node, bool isrooted, map<unsigned long long, unsigned long long> &bipcount, int depth)
{
    depth++;
    if( (isrooted && depth > 1 && node->Nchildren > 0) ||
        (! isrooted && depth > 2 && node->Nchildren > 0))
    {
        map<unsigned long long, unsigned long long>::iterator iter;
        iter = bipcount.find(node->hv2);
        if(iter != bipcount.end())
        {
            bipcount[node->hv2]++;
        } else
        {
            bipcount[node->hv2] = 1;
        }
    }

    for(int i = 0; i < node->Nchildren; i++)
    {
        bipartcount(node->child[i], isrooted, bipcount, depth);
    }
}

bool TreeOPE::buildconsensustree(NEWICKTREE *&tree,
                                 const Array<double> &confreq,
                                 const Array<unsigned long long> &conhash,
                                 Array<char> *contreebitstr,
                                 unsigned int conbipnum,
                                 int bitlength)
{
    // check compatible
    bool firstcontainsecond;
    Array<char> *flipbitstr = new Array<char> [conbipnum];
    for(int i = 0; i < conbipnum; i++)
    {
        flipbitstr[i] = contreebitstr[i];
        flipbitstr[i].Flipbit(bitlength);
    }

    for(int i = 0; i < conbipnum; i++)
    {
        for(int j = 0; j < i; j++)
        {
            if(  !Array<char>::bitAcontainB(contreebitstr[i], contreebitstr[j], bitlength, firstcontainsecond)
              && !Array<char>::bitAcontainB(contreebitstr[i], flipbitstr[j], bitlength, firstcontainsecond)
              && !Array<char>::bitAcontainB(flipbitstr[i], contreebitstr[j], bitlength, firstcontainsecond)
              && !Array<char>::bitAcontainB(flipbitstr[i], flipbitstr[j], bitlength, firstcontainsecond)
                 )
            {
                std::cout << i << "-th tree and " << j << "-th tree are not compatible!\n\n";
                /*
                contreebitstr[i].printbits(bitlength);
                cout << endl;
                flipbitstr[i].printbits(bitlength);
                cout << endl;
                contreebitstr[j].printbits(bitlength);
                cout << endl;
                flipbitstr[j].printbits(bitlength);
                cout << endl;
                if(Array<char>::bitAcontainB(flipbitstr[i], contreebitstr[j], bitlength, firstcontainsecond))
                    cout << "true" << endl;
                else
                    cout << "false" << endl;//--
                */
                return false;
            }
        }
    }
    delete [] flipbitstr;


    // build consensus tree
    tree = new NEWICKTREE;
    tree->root = new NEWICKNODE;
    tree->root->bitstr = NULL;
    tree->root->hv1 = 0;
    tree->root->hv2 = 0;
    tree->root->label = NULL;
    tree->root->parent = NULL;
    tree->root->weight = 0;
    tree->root->Nchildren = bitlength;
    tree->root->child = new NEWICKNODE *[bitlength];

    int capacity = 32; //consistent with readplainlabel in TreeOPE.cpp
    for(int i = 0; i < bitlength; i++)
    {
        tree->root->child[i] = new NEWICKNODE;
        tree->root->child[i]->bitstr = NULL;
        tree->root->child[i]->child = NULL;
        tree->root->child[i]->hv1 = 0;
        tree->root->child[i]->hv2 = 0;
        stringstream ss;
        string stdstr;
        ss << i + 1;
        stdstr = ss.str();
        tree->root->child[i]->label = new char [capacity];
        strcpy(tree->root->child[i]->label, stdstr.c_str());
        tree->root->child[i]->Nchildren = 0;
        tree->root->child[i]->parent = tree->root;
        tree->root->child[i]->weight = 1;
    }

    for(int i = 0; i < conbipnum; i++)//1; i++)//--
    {
        bool ischecked;
        Addbipart(tree->root, confreq[i], conhash[i], contreebitstr[i], bitlength, ischecked);
    }
    return true;
}

bool TreeOPE::Addbipart(NEWICKNODE* startNode, double freq, unsigned long long hash, Array<char> &bitstr, int NumTaxa, bool &iscontained)
{
    int Bytelength = (int) ((NumTaxa + 7) / 8);

    if (startNode->Nchildren == 0)
    {
        unsigned int idx = atoi(startNode->label) - 1;

        if(startNode->bitstr != NULL)
            delete startNode->bitstr;

        startNode->bitstr = new Array<char> (Bytelength);

        if(! startNode->bitstr->SetBitArray(idx + 1))
            std::cout << "Warning: Errors in setting bit string!\n\n";

        Array<char>::bitAcontainB(bitstr, *(startNode->bitstr), NumTaxa, iscontained);
/*
        bitstr.printbits(NumTaxa);
        cout << endl;
        startNode->bitstr->printbits(NumTaxa);
        cout << endl;//--
        std::cout << "iscontained:" << iscontained << std::endl;//---
*/
        return true;
    } else
    {
        bool containrelation;
        int numtrues = 0;
        int numflase = 0;
        int *trueidx = new int [startNode->Nchildren];
        int *flaseidx = new int [startNode->Nchildren];
        for (int i = 0; i < startNode->Nchildren; i++)
        {
            if(! Addbipart(startNode->child[i], freq, hash, bitstr, NumTaxa, containrelation))
            {
                delete [] trueidx;
                delete [] flaseidx;
                return false;
            }
            if(containrelation)
            {
                trueidx[numtrues] = i;
                numtrues++;
            } else
            {
                flaseidx[numflase] = i;
                numflase++;
            }
        }
//        std::cout << "numtrues:" << numtrues << ", numflase:" << numflase << std::endl;//---

        if(startNode->bitstr != NULL)
            delete startNode->bitstr;
        startNode->bitstr = new Array<char> (Bytelength);

        for (int i = 0; i < startNode->Nchildren; i++)
        {
            startNode->bitstr->ORbitOPE(*(startNode->bitstr), *(startNode->child[i]->bitstr));
        }

        Array<char>::bitAcontainB(bitstr, *(startNode->bitstr), NumTaxa, iscontained);
        if(numtrues == 0 || numtrues == startNode->Nchildren)
        {
            delete [] trueidx;
            delete [] flaseidx;
            return true;
        }

        NEWICKNODE *node = new NEWICKNODE;
        node->bitstr = NULL;
        node->hv1 = 0;
        node->hv2 = hash;
        node->label = NULL;
        node->weight = freq;
        node->Nchildren = numtrues;
        node->child = new NEWICKNODE *[numtrues];
        node->parent = startNode;
        for(int i = 0; i < numtrues; i++)
        {
            node->child[i] = startNode->child[trueidx[i]];
        }

        for(int i = 0; i < numflase; i++)
        {
            startNode->child[i] = startNode->child[flaseidx[i]];
        }
        startNode->child[numflase] = node;
        startNode->Nchildren = numflase + 1;

        delete [] trueidx;
        delete [] flaseidx;
        return false;
    }
}

#endif
