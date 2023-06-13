#pragma once

#include <iostream>

#include <stack>
#include <queue>
#include <cmath>

#include "lexer.hpp"
#include "nwq_util.hpp"

using namespace std;
using namespace lexertk;

int get_precedence(token::token_type op)
{
    switch (op)
    {
    case token::e_add:
    case token::e_sub:
        return 0;

    case token::e_mul:
    case token::e_div:
        return 1;
    case token::e_pow:
        return 2;
    case token::e_negative:
        return 3;

    default:
        return -1;
    }
}

int compare_operators(token::token_type op1, token::token_type op2)
{
    int op1_precedence = get_precedence(op1);
    int op2_precedence = get_precedence(op2);

    if (op1_precedence == -1 || op2_precedence == -1)
        throw runtime_error("Error: Unrecognized operator");

    return op1_precedence > op2_precedence ? 1 : op1_precedence < op2_precedence ? -1
                                                                                 : 0;
}

bool is_operator(token::token_type op)
{
    return (op == token::e_pow) ||
           (op == token::e_mul) ||
           (op == token::e_div) ||
           (op == token::e_add) ||
           (op == token::e_sub) ||
           (op == token::e_negative);
}

/**
 * Parse the expression token using Shunting Yard Algorithm
 */
double parse_expr(vector<token> tokens, int start, int end)
{
    stack<token> op_stack;
    queue<token> op_queue;

    for (int i = start; i < end; i++)
    {
        token t = tokens[i];

        switch (t.type)
        {
        case token::e_number:
        case token::e_pi:
            op_queue.push(t);
            break;

        case token::e_func:
            op_stack.push(t);
            break;

        case token::e_pow:

        case token::e_lbracket:
            op_stack.push(t);
            break;

        case token::e_rbracket:
            while (op_stack.top().type != token::e_lbracket)
            {
                op_queue.push(op_stack.top());
                op_stack.pop();
            }
            op_stack.pop();

            if (!op_stack.empty() && op_stack.top().type == token::e_func)
            {
                op_queue.push(op_stack.top());
                op_stack.pop();
            }
            break;

        case token::e_sub:
            if (i == 0 || !(tokens[i - 1].type == token::e_number || tokens[i - 1].type == token::e_pi || tokens[i - 1].type == token::e_rbracket))
            {
                t.type = token::e_negative;
                t.value = "NEG";
                op_stack.push(t);
                break;
            }
        case token::e_add:
        case token::e_mul:
        case token::e_div:
            while (!op_stack.empty() && is_operator(op_stack.top().type))
            {
                if (compare_operators(t.type, op_stack.top().type) <= 0)
                {
                    op_queue.push(op_stack.top());
                    op_stack.pop();
                }
                else
                {
                    break;
                }
            }

            op_stack.push(t);
            break;

        default:
            cout << "Unknown token type: " << t.to_str(t.type).c_str() << " " << t.value << endl;
            break;
        }
    }

    while (!op_stack.empty())
    {
        op_queue.push(op_stack.top());
        op_stack.pop();
    }

    // int num_token = op_queue.size();
    // while (num_token > 0)
    // {
    //     auto t = op_queue.front();

    //     cout << t.value << " , ";
    //     op_queue.pop();
    //     op_queue.push(t);
    //     num_token--;
    // }
    // cout << endl;

    stack<double> val_stack;
    double val1, val2;
    while (!op_queue.empty())
    {
        token t = op_queue.front();
        op_queue.pop();

        switch (t.type)
        {
        case token::e_number:
            val_stack.push(stod(t.value));
            // cout << "Pushed " << t.value << endl;
            break;

        case token::e_pi:
            val_stack.push(PI);
            // cout << "Pushed PI" << endl;
            break;

        case token::e_add:
            val2 = val_stack.top();
            val_stack.pop();

            if (!val_stack.empty())
            {
                val1 = val_stack.top();
                val_stack.pop();
            }
            else
            {
                val1 = 0;
            }

            val_stack.push(val1 + val2);

            // cout << val1 << " + " << val2 << endl;
            break;

        case token::e_sub:
            val2 = val_stack.top();
            val_stack.pop();

            if (val_stack.empty())
                val_stack.push(-val2);
            else
            {
                val1 = val_stack.top();

                val_stack.pop();
                val_stack.push(val1 - val2);
            }
            // cout << val1 << " - " << val2 << endl;
            break;

        case token::e_mul:
            val2 = val_stack.top();
            val_stack.pop();
            val1 = val_stack.top();
            val_stack.pop();
            val_stack.push(val1 * val2);

            // cout << val1 << " * " << val2 << endl;
            break;

        case token::e_div:
            val2 = val_stack.top();
            val_stack.pop();
            val1 = val_stack.top();
            val_stack.pop();
            val_stack.push(val1 / val2);

            // cout << val1 << " / " << val2 << endl;
            break;

        case token::e_pow:
            val2 = val_stack.top();
            val_stack.pop();
            val1 = val_stack.top();
            val_stack.pop();
            val_stack.push(pow(val1, val2));

            // cout << val1 << " ^ " << val2 << endl;
            break;

        case token::e_func:
        case token::e_negative:

            val1 = val_stack.top();
            val_stack.pop();

            if (t.value == "cos")
            {
                val_stack.push(cos(val1));
                // cout << "cos " << val1 << endl;
            }
            else if (t.value == "sin")
            {
                val_stack.push(sin(val1));
                // cout << "sin " << val1 << endl;
            }
            else if (t.type == token::e_negative)
            {
                val_stack.push(-val1);
            }
            else
                cout << "UNRECOGNIZED FUNCTION: " << t.value << endl;
            break;

        default:
            cout << "UNRECOGNIZED TOKEN: " << t.value << endl;
            break;
        }
    }

    return val_stack.top();
}