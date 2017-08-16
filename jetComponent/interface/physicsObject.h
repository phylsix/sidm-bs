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
    
    class Jet : public objBase{
    };

    class Ep : public objBase{
    };
}

#endif
