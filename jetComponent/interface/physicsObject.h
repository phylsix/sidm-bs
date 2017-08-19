#ifndef SIDMBS_PHYSICSOBJECT_H
#define SIDMBS_PHYSICSOBJECT_H

namespace sidm {

    class objBase{
    public:
        int _eventId;
        float _pt;
        float _eta;
        float _phi;
        float _energy;
        float _mass;
    };

    class objCompound : public objBase{
    public:
        float _invM;
        float _dEta;
        float _dPhi;
        float _dR;

    };
    
    class Jet : public objBase{
    };

    class Ep : public objBase{
    };

    class Zp : public objCompound{
    public:
        Ep p;
        Ep e;
        float _dv_x;
        float _dv_y;
        float _dv_z;
    };

    class Ps : public objCompound{
    public:
        Zp zp_0;
        Zp zp_1;
    };
}

#endif
