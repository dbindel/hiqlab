/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef MODEL_CIRCUIT_H
#define MODEL_CIRCUIT_H

#include "element.h"
#include "qmatrix.h"


/** Scalar Resistor element.
 */
class Resistor : public Element {
 public:

    /** Construct a new resistor element
     *
     * @param resist  Resistance
     */
    Resistor(double resist, double resist_i = 0.0);

    ~Resistor();

    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double resist, resist_i;
};


/** Scalar Capacitor element.
 */
class Capacitor : public Element {
 public:

    /** Construct a new capacitor element
     *
     * @param cap  Capacitance
     */
    Capacitor(double cap, double cap_i = 0.0);

    ~Capacitor();

    void initialize(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double cap, cap_i;
};


/** Scalar Capacitor element. with ONE Lagrange multiplier (charge)
 */
class Capacitor1 : public Element {
 public:

    /** Construct a new capacitor element
     *
     * @param cap  Capacitance
     */
    Capacitor1(double cap, double cap_i = 0.0);

    ~Capacitor1();

    void initialize(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double cap, cap_i;
};


/** Scalar Capacitor element. With TWO LAGRANGE MULTIPLIERS(charge, current)
 */
class Capacitor2 : public Element {
 public:

    /** Construct a new capacitor element
     *
     * @param cap  Capacitance
     */
    Capacitor2(double cap, double cap_i = 0.0);

    ~Capacitor2();

    void initialize(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double cap, cap_i;
};


/** Scalar Inductor element.
 */
class Inductor : public Element {
 public:

    /** Construct a new inductor element
     *
     * @param induct  Inductance
     */
    Inductor(double induct, double induct_i = 0.0);

    ~Inductor();

    void initialize(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    double induct, induct_i;
};


/** Voltage or current source element.
 */
class VIsrc : public Element {
 public:

    /** Construct a new Voltage or current source element
     *
     */
    VIsrc();

    ~VIsrc();

    void initialize(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

};


/** Electrode element.
 */
class Electrode : public Element {
 public:

    /** Construct a new Electrode element
     *
     * @param L     Lua state
     * @param func  Lua function for global variable
     * @param enode Node to connect the electrode
     */

    Electrode(int vglobalid, double lt=1);

    ~Electrode();

    void initialize(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    int  vglobalid;
    double lt;

};


/** Electrode element. TWO LAGRANGE MULTIPLIERS
 */
class Electrode2 : public Element {
 public:

    /** Construct a new Electrode element
     *
     * @param L     Lua state
     * @param func  Lua function for global variable
     * @param enode Node to connect the electrode
     */

    Electrode2(int vglobalid, double lt=1);

    ~Electrode2();

    void initialize(Mesh* mesh, int eltid);
    void assemble_dR(Mesh* mesh, int eltid, QAssembler* K,
                     double cx=1, double cv=0, double ca=0);
    void assemble_R(Mesh* mesh, int eltid);

 private:
    int  vglobalid;
    double lt;

};


#endif /* MODEL_CIRCUIT_H */
