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

    };
    
    class Jet : public objBase{
    };

    class Ep : public objBase{
    };

    class Zp : public objCompound{
    public:
        Ep p;
        Ep e;
        float _dR_ep;
    };

    class Ps : public objCompound{
    public:
        Zp zp_0;
        Zp zp_1;
    };
}

#endif
