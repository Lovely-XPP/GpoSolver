/*
    GpoSolver, Library for Global Polynomial Optimization
    Copyright (C) 2014-2015 Jan Heller, <hellej1@cmp.felk.cvut.cz>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GPOPROBLEM_H
#define GPOPROBLEM_H

#include <iostream>
#include <set>
#include <map>
#include <list>

#include <Eigen/Dense>

namespace GpoSolver
{
    using namespace std;

    class GpoProblem
    {
    public:
        // SDP objective non-zero vector entry
        typedef struct
        {
            unsigned int mon;
            double val;
        } objData;

        // SDP objective variable vector entry
        typedef struct
        {
            unsigned int mon;
        } objVar;

        // SDP constraints non-zero matrix entry
        typedef struct
        {
            unsigned int con;
            unsigned int mon;
            unsigned int i;
            unsigned int j;
            double val;
        } conData;

        // SDP constraints matrix variable
        typedef struct
        {
            int con;
            int mon;
            int i;
        } conVar;

        // SDP constraints matrix block info
        typedef struct
        {
            unsigned int con;
            unsigned int mon;
            unsigned int num_entries;
        } conBlock;

        // Problem status
        typedef enum
        {
            OK = 0,
        } Status;

        // Virtual methods to be automatically generated ==========================

        virtual int getNumBlocks(void) const = 0;
        virtual int getNumConstraints(void) const = 0;
        virtual int getNumProblemVariables(void) const = 0;
        virtual int getNumPolyConstraints(void) const = 0;
        virtual int getNumObjectiveParams(void) const = 0;
        virtual int getNumConstraintsParams(void) const = 0;
        virtual int getProblemSize(void) const = 0;
        virtual int getRelaxationOrder(void) const = 0;
        virtual int getMaxDegree(void) const = 0;
        virtual int getRankShift(void) const = 0;
        virtual int getBlockSize(int) const = 0;
        virtual const int *getBlockSizes(void) const = 0;

        virtual double evalPolyObjective(const double *, const double *, const double *) const = 0;
        virtual void evalPolyConstraints(const double *, const double *, double *) const = 0;

        virtual const objData *getSdpObjectiveData(void) const = 0;
        virtual const objVar *getSdpObjectiveVariables(void) const = 0;
        virtual int getNumSdpObjectiveData(void) const = 0;
        virtual int getNumSdpObjectiveVariables(void) const = 0;
        virtual void evalSdpObjectiveVariables(const double *, const double *, double **) const = 0;

        virtual const conData *getSdpConstraintsData(void) const = 0;
        virtual const conVar *getSdpConstraintsVariables(void) const = 0;
        virtual int getNumSdpConstraintsData(void) const = 0;
        virtual int getNumSdpConstraintsVariables(void) const = 0;
        virtual void evalSdpConstraintsVariables(const double *, const double *, double **) const = 0;

        // SDP problem construction ===============================================

        typedef map<pair<int, int>, pair<const conData *, int>> blockMap;
        typedef map<int, set<int>> conBlockMap;
        typedef map<pair<int, int>, list<pair<int, int>>> conVarMap;

        blockMap block_map;
        conBlockMap conblock_map;
        conVarMap convars_map;

        conData _con_data[10000] = {};
        conVar _con_vars[10000] = {};
        objData _obj_data[1000] = {};
        objVar _obj_vars[1000] = {};
        int _block_sizes[100] = {};
        int _num_blocks = 3;
        int _num_cons = 9;
        int _num_pcons = 2;
        int _num_pvars = 3;
        int _num_rpars = 0;
        int _num_ppars = 0;
        int _problem_size = 0;
        int _relax_order = 1;
        int _max_degree = 2;
        int _rank_shift = 1;
        int _num_obj_data = 1;
        int _num_obj_vars = 0;
        int _num_con_data = 0;
        int _num_con_vars = 0;

        // user define data
        int mat_size = 0;

        int getNumNonzeroBlocks(void) const
        {
            return block_map.size();
        }

        void getBlock(int i, int j, const conData *&data, int &ndata) const
        {
            blockMap::const_iterator blk = block_map.find(make_pair(i, j));
            if (blk == block_map.end())
            {
                data = NULL;
                ndata = 0;
            }
            else
            {
                data = blk->second.first;
                ndata = blk->second.second;
            }
        }

        const conData *getBlockData(int i, int j) const
        {
            int ndata;
            const GpoProblem::conData *data;
            getBlock(i, j, data, ndata);
            return data;
        }

        int getNumBlockEntries(int i, int j) const
        {
            int ndata;
            const conData *data;
            getBlock(i, j, data, ndata);
            return ndata;
        }

        const set<int> *getBlocksByConstraint(int i) const
        {
            conBlockMap::const_iterator con = conblock_map.find(i);
            if (con == conblock_map.end())
                return NULL;
            else
                return &(con->second);
        }

        const list<pair<int, int>> *getVarsByBlock(int i, int j) const
        {
            conVarMap::const_iterator blk = convars_map.find(make_pair(i, j));
            if (blk == convars_map.end())
                return NULL;
            else
                return &(blk->second);
        }

        virtual ~GpoProblem() {}

        void init(void)
        {
            _num_blocks = getNumBlocks();
            _num_cons = getNumConstraints();
            _num_pcons = getNumPolyConstraints();
            _num_pvars = getNumProblemVariables();
            _num_rpars = getNumObjectiveParams();
            _num_ppars = getNumConstraintsParams();
            _problem_size = getProblemSize();
            _relax_order = getRelaxationOrder();
            _max_degree = getMaxDegree();
            _rank_shift = getRankShift();
            _num_obj_data = getNumSdpObjectiveData();
            _num_obj_vars = getNumSdpObjectiveVariables();
            _num_con_data = getNumSdpConstraintsData();
            _num_con_vars = getNumSdpConstraintsVariables();

            const int *block_size = getBlockSizes();
            for (size_t i = 0; i < _num_blocks; i++)
            {
                _block_sizes[i] = block_size[i];
            }

            const conData *con_data = getSdpConstraintsData();
            for (int i = 0; i < _num_con_data; i++)
            {
                _con_data[i] = con_data[i];
            }

            const conVar *con_vars = getSdpConstraintsVariables();
            for (int i = 0; i < _num_con_vars; i++)
            {
                _con_vars[i] = con_vars[i];
            }

            const objData *obj_data = getSdpObjectiveData();
            for (int i = 0; i < _num_obj_data; i++)
            {
                _obj_data[i] = obj_data[i];
            }

            const objVar *obj_vars = getSdpObjectiveVariables();
            for (int i = 0; i < _num_obj_vars; i++)
            {
                _obj_vars[i] = obj_vars[i];
            }

            const conData *data = _con_data;
            block_map.clear();
            for (int i = 0; i < _num_con_data; i++)
            {
                pair<int, int> key = make_pair(data[i].con, data[i].mon);
                blockMap::iterator blk = block_map.find(key);
                if (blk == block_map.end())
                    block_map[key] = make_pair(data + i, 1);
                else
                    blk->second.second++;

                conBlockMap::iterator mon = conblock_map.find(data[i].mon);
                if (mon == conblock_map.end())
                {
                    conblock_map[data[i].mon] = set<int>();
                    conblock_map[data[i].mon].insert(data[i].con);
                }
                else
                    mon->second.insert(data[i].con);
            }

            const GpoProblem::conVar *vars = _con_vars;
            convars_map.clear();
            for (int i = 0; i < _num_con_vars; i++)
            {
                pair<int, int> key = make_pair(vars[i].con, vars[i].mon);
                pair<int, int> val = make_pair(i, vars[i].i);
                conVarMap::iterator blk = convars_map.find(key);
                if (blk == convars_map.end())
                {
                    convars_map[key] = list<pair<int, int>>();
                    convars_map[key].push_back(val);
                }
                else
                    blk->second.push_back(val);
            }
        }

        void init(const int n)
        {
            mat_size = n;
            _num_blocks = 3;
            _num_cons = 9;
            _num_pcons = 2;
            _num_pvars = 3;
            _num_rpars = 0;
            _num_ppars = 0;
            _problem_size = 5 + n;
            _relax_order = 1;
            _max_degree = 2;
            _rank_shift = 1;
            _num_obj_data = 1;
            _num_obj_vars = 0;
            _num_con_data = 13 + n * (n + 1) + n;
            _num_con_vars = n * (n + 1);

            const int *block_size = getBlockSizes();
            for (size_t i = 0; i < _num_blocks; i++)
            {
                _block_sizes[i] = block_size[i];
            }
            _block_sizes[2] = n;

            const conData *con_data = getSdpConstraintsData();
            int const_size = 13;
            int m = const_size;
            for (int i = 0; i < const_size; i++)
            {
                _con_data[i] = con_data[i];
            }
            for (unsigned int i = 0; i < n; i++)
            {
                for (unsigned int j = i; j < n; j++)
                {
                    _con_data[m] = {2, 1, i, j, 0.0};
                    m++;
                }
            }
            for (unsigned int i = 0; i < n; i++)
            {
                for (unsigned int j = i; j < n; j++)
                {
                    _con_data[m] = {2, 2, i, j, 0.0};
                    m++;
                }
            }
            for (unsigned int i = 0; i < n; i++)
            {
                _con_data[m] = {2, 3, i, i, 1.0};
                m++;
            }

            const conVar *con_vars = getSdpConstraintsVariables();
            for (int i = 0; i < (n * (n + 1) / 2); i++)
            {
                _con_vars[i] = {2, 1, i};
            }
            for (int i = 0; i < (n * (n + 1) / 2); i++)
            {
                _con_vars[i + n * (n + 1) / 2] = {2, 2, i};
            }

            const objData *obj_data = getSdpObjectiveData();
            for (int i = 0; i < _num_obj_data; i++)
            {
                _obj_data[i] = obj_data[i];
            }

            const objVar *obj_vars = getSdpObjectiveVariables();
            for (int i = 0; i < _num_obj_vars; i++)
            {
                _obj_vars[i] = obj_vars[i];
            }

            const conData *data = _con_data;
            block_map.clear();
            for (int i = 0; i < _num_con_data; i++)
            {
                pair<int, int> key = make_pair(data[i].con, data[i].mon);
                blockMap::iterator blk = block_map.find(key);
                if (blk == block_map.end())
                    block_map[key] = make_pair(data + i, 1);
                else
                    blk->second.second++;

                conBlockMap::iterator mon = conblock_map.find(data[i].mon);
                if (mon == conblock_map.end())
                {
                    conblock_map[data[i].mon] = set<int>();
                    conblock_map[data[i].mon].insert(data[i].con);
                }
                else
                    mon->second.insert(data[i].con);
            }

            const GpoProblem::conVar *vars = _con_vars;
            convars_map.clear();
            for (int i = 0; i < _num_con_vars; i++)
            {
                pair<int, int> key = make_pair(vars[i].con, vars[i].mon);
                pair<int, int> val = make_pair(i, vars[i].i);
                conVarMap::iterator blk = convars_map.find(key);
                if (blk == convars_map.end())
                {
                    convars_map[key] = list<pair<int, int>>();
                    convars_map[key].push_back(val);
                }
                else
                    blk->second.push_back(val);
            }
        }
    }; // GpoProblem

} // GpoSolver namespace

#endif // GPOPROBLEM
