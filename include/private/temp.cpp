inline void registerGates()
{

    GateFactory::getInstance().registerGate(OP::X, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {0, 1, 1, 0};
                                                 ValType imag[4] = {0, 0, 0, 0};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::Y, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {0, 0, 0, 0};
                                                 ValType imag[4] = {0, -1, 1, 0};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::Z, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {1, 0, 0, -1};
                                                 ValType imag[4] = {0, 0, 0, 0};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::H, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {S2I, S2I, S2I, -S2I};
                                                 ValType imag[4] = {0, 0, 0, 0};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::S, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {1, 0, 0, 0};
                                                 ValType imag[4] = {0, 0, 0, 1};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::SDG, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {1, 0, 0, 0};
                                                 ValType imag[4] = {0, 0, 0, -1};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::T, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {1, 0, 0, S2I};
                                                 ValType imag[4] = {0, 0, 0, S2I};
                                                  memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::TDG, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {1, 0, 0, S2I};
                                                 ValType imag[4] = {0, 0, 0, -S2I};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::RI, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                 ValType real[4] = {cos(theta), 0, 0, cos(theta)};
                                                 ValType imag[4] = {sin(theta), 0, 0, sin(theta)};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::RX, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                 ValType real[4] = {cos(HALF * theta), 0, 0, cos(HALF * theta)};
                                                 ValType imag[4] = {0, -sin(HALF * theta), -sin(HALF * theta), 0};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::RY, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                 ValType real[4] = {cos(HALF * theta), -sin(HALF * theta), sin(HALF * theta), cos(HALF * theta)};
                                                 ValType imag[4] = {0};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::RZ, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                 ValType real[4] = {cos(HALF * theta), 0, 0, cos(HALF * theta)};
                                                 ValType imag[4] = {-sin(HALF * theta), 0, 0, sin(HALF * theta)};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::SX, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[4] = {HALF, HALF, HALF, HALF};
                                                 ValType imag[4] = {HALF, -HALF, -HALF, HALF};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::P, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                 ValType real[4] = {1, 0, 0, cos(theta)};
                                                 ValType imag[4] = {0, 0, 0, sin(theta)};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::U, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType phi, ValType lam)
                                            {
                                                 ValType real[4] = {cos(HALF * theta),
                                                                    -cos(lam) * sin(HALF * theta),
                                                                    cos(phi) * sin(HALF * theta),
                                                                    cos(phi + lam) * cos(HALF * theta)};
                                                 ValType imag[4] = {0,
                                                                    -sin(lam) * sin(HALF * theta),
                                                                    sin(phi) * sin(HALF * theta),
                                                                    sin(lam + phi) * cos(HALF * theta)};
                                                 memcpy(gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 4 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::CX, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 0, 1,
                                                                     0, 0, 1, 0};
                                                 ValType imag[16] = {0};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::CY, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, -1,
                                                                     0, 0, 1, 0};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::CZ, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, -1};
                                                 ValType imag[16] = {0};
                                                   memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::CH, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, S2I, S2I,
                                                                     0, 0, S2I, -S2I};
                                                 ValType imag[16] = {0};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::CS, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, 0};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 1};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::CSDG, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, 0};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, -1};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::CT, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, S2I};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, S2I};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::CTDG, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                                ValType real[16] = {1, 0, 0, 0,
                                                0, 1, 0, 0,
                                                0, 0, 1, 0,
                                                0, 0, 0, S2I};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                0, 0, 0, 0, 
                                                0, 0, 0, 0,
                                                0, 0, 0, -S2I};  
                                                memcpy(gm_real, real, 16 * sizeof(ValType));
                                                memcpy(gm_imag, imag, 16 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::CRX, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, cos(HALF * theta), 0,
                                                                     0, 0, 0, cos(HALF * theta)};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, -sin(HALF * theta),
                                                                     0, 0, -sin(HALF * theta), 0};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::CRY, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, cos(HALF * theta), -sin(HALF * theta),
                                                                     0, 0, sin(HALF * theta), cos(HALF * theta)};
                                                 ValType imag[16] = {0};
                                                 memcpy(gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(gm_imag, imag, 16 * sizeof(ValType)); });
    GateFactory::getInstance().registerGate(OP::CRZ, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                ValType real[16] = {1, 0, 0, 0,
                                                                    0, 1, 0, 0,
                                                                    0, 0, cos(HALF * theta), 0,
                                                                    0, 0, 0, cos(HALF * theta)};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, -sin(HALF * theta), 0,
                                                                    0, 0, 0, sin(HALF * theta)};
                                                    memcpy(gm_real, real, 16 * sizeof(ValType));
                                                    memcpy(gm_imag, imag, 16 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::CSX, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType, ValType, ValType)
                                            {
                                               ValType real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, HALF, HALF,
                                   0, 0, HALF, HALF};
            ValType imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, HALF, -HALF,
                                   0, 0, -HALF, HALF};
                                                    memcpy(gm_real, real, 16 * sizeof(ValType));
                                                    memcpy(gm_imag, imag, 16 * sizeof(ValType)); });

    GateFactory::getInstance().registerGate(OP::CP, [](ValType gm_real[], ValType gm_imag[], SIM_TYPE, ValType theta, ValType, ValType)
                                            {
                                                ValType real[16] = {1, 0, 0, 0,
                           0, 1, 0, 0,
                           0, 0, 1, 0,
                           0, 0, 0, cos(theta)};
    ValType imag[16] = {0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, 0,
                           0, 0, 0, sin(theta)};
                                                    memcpy(gm_real, real, 16 * sizeof(ValType));
                                                    memcpy(gm_imag, imag, 16 * sizeof(ValType)); });
}