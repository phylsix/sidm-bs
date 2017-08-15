#ifndef SIDMBS_PHYSICSOBJECT_H
#define SIDMBS_PHYSICSOBJECT_H

namespace sidm {
    
    struct Jet{
        int _eventId;
        float _pt;
        float _eta;
        float _phi;
        float _energy;
        float _mass;
    };
}

#endif
