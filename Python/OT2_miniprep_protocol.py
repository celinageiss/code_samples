# =============================================================================
# Zyppy Miniprep protocol
# =============================================================================

# Set your parameters
start_column = 1  # Start column on plate (value between 1 and 6)
total_columns = 1  # Put number of columns here (value between 1 and 6)
drying_time = 30  # min
elution_buffer_volume = 50  # µl
elution_column = 1

# =============================================================================
# Now load the protocol on the robot. Please do NOT change the script below!
# =============================================================================

# imports
from opentrons import protocol_api, types

# metadata
metadata = {
    'protocolName': 'Plasmid DNA Purification',
    'author': 'Celina Geiss <celina.geiss@dkfz-heidelberg.de>',
    'description': 'Purification of plasmid DNA using the Zymo Zyppy-96 Plasmid MagBead Miniprep Kit.',
    'apiLevel': '2.6'
}


################# Labware #################
# 10. Tiprack   11. Tiprack  12. TRASH
#  7. Tiprack    8. Tiprack   9.
#  4. Waste      5. Tiprack   6.
#  1. Magnet     2. Reservoir 3. Elution
###########################################

# ----------------
# Define functions
# ----------------

# Pick up tip if tip is not yet attached
def check_tip_attached(pipette):
    if not pipette.hw_pipette['has_tip']:
        pipette.pick_up_tip()
    return ()


# Mix from top to bottom and back
def mix_up_down(pipette,
                vol,
                plate,
                col,
                top_height,
                bottom_height,
                steps,  # steps in mm
                speed_asp=300,
                speed_disp=600):

    pipette.flow_rate.aspirate = speed_asp
    pipette.flow_rate.dispense = speed_disp
    # mix from top to bottom
    [pipette.mix(1,
                 vol,
                 plate[col][0].bottom(z=k)) for k in reversed(range(bottom_height, top_height, steps))]
    # mix from bottom to top
    [pipette.mix(1,
                 vol,
                 plate[col][0].bottom(z=k)) for k in range(bottom_height, top_height, steps)]
    pipette.blow_out(plate[col][0])
    pipette.touch_tip(plate[col][0], v_offset=-2)


# Transfer larger volumes
def transfer_large_vol(pipette,
                       total_vol,
                       source_well,
                       dest_plate,
                       col,
                       speed_asp=150,
                       tip_vol=200):

    pipette.flow_rate.aspirate = speed_asp
    vol_remaining = total_vol
    while vol_remaining > 0:
        transfer_vol = min(vol_remaining, (tip_vol-10))
        pipette.aspirate(transfer_vol,
                         source_well)
        pipette.air_gap(10)
        pipette.dispense(transfer_vol + 10, dest_plate[col][0])
        pipette.blow_out(dest_plate[col][0].top(z=-3))
        vol_remaining -= transfer_vol
    pipette.touch_tip(dest_plate[col][0], v_offset=-2)


# Define x-offset for odd or even columns
def odd_or_even(plate,
                col,
                x_offset=0.75,
                z_offset=0):
    side = col % 2
    offset = -((-1) ** side) * x_offset
    center_loc = plate[col][0].bottom()
    offset_loc = center_loc.move(types.Point(x=offset, y=0, z=z_offset))
    return (offset_loc)


# Transfer supernatant to waste
def discard_supernatant(pipette,
                       total_vol,
                       plate,
                       col,
                       waste,
                       bottom_offset,
                       tip_vol=200,
                       speed_asp=30):

    pipette.flow_rate.aspirate = speed_asp
    vol_remaining = total_vol
    pipette.pick_up_tip()
    while vol_remaining > 0:
        transfer_vol = min(vol_remaining, (tip_vol-10))
        if vol_remaining <= 190:
            plate_offset = odd_or_even(plate, col, z_offset=bottom_offset-1.25)  # determine offset
        else:
            plate_offset = odd_or_even(plate, col, z_offset=2)  # determine offset
        pipette.aspirate(transfer_vol,
                         plate_offset)
        pipette.air_gap(10)
        pipette.dispense(transfer_vol + 10, waste)
        pipette.blow_out()
        vol_remaining -= transfer_vol
    pipette.drop_tip()

# -----------------
# Protocol function
# -----------------

def run(protocol: protocol_api.ProtocolContext):
    # --------------
    # Define labware
    # --------------

    magmol = protocol.load_module('Magnetic Module', 1)
    deep_plate = magmol.load_labware('usascientific_96_wellplate_2.4ml_deep').columns()
    reservoir = protocol.load_labware('usascientific_12_reservoir_22ml', 2)
    elution_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 3).columns()
    waste = protocol.load_labware('agilent_1_reservoir_290ml', 4).wells()[0]

    # Tipracks
    tipslots = [5, 7, 8, 10, 11]
    tipracks = [protocol.load_labware('opentrons_96_tiprack_300ul', slot) for slot in tipslots]

    # Pipettes
    #     p300_single = protocol.load_instrument('p300_single','left', tip_racks=tipracks)
    p300_multi = protocol.load_instrument('p300_multi', 'right', tip_racks=tipracks)

    # Buffers in 12-well reservoir
    lysis_buffer = reservoir.wells()[0]
    neutralisation_1 = reservoir.wells()[1]
    neutralisation_2 = reservoir.wells()[2]
    neutralisation = [neutralisation_1, neutralisation_2]
    clearing_beads = reservoir.wells()[3]
    binding_beads = reservoir.wells()[4]
    endo_wash = reservoir.wells()[5]
    zyppy_wash_1 = reservoir.wells()[6]
    zyppy_wash_2 = reservoir.wells()[7]
    zyppy_wash_3 = reservoir.wells()[8]
    zyppy_wash = [zyppy_wash_1, zyppy_wash_2, zyppy_wash_3]
    elution_buffer = reservoir.wells()[9]

    # Define sample wells
    start_col = start_column - 1
    end_col = start_col + total_columns
    cols = range(start_col, end_col)
    elution_col = elution_column - 1

    # Magnet engage height
    mag_engage_height = 13
    magmol.disengage()

    # -----------------
    # Start of protocol
    # -----------------

    # TODO: Return tips after mixing

    # Set parameters for aspiration and dispense
    p300_multi.well_bottom_clearance.dispense = 39  # dispense 1 mm from top of deep_plate

    # Transfer Deep-Blue Lysis Buffer
    protocol.comment("Dispensing Blue Lysis Buffer.")
    for col in cols:
        p300_multi.pick_up_tip()
        p300_multi.transfer(100, lysis_buffer, deep_plate[col][0], new_tip='never')
        mix_up_down(pipette=p300_multi,
                    vol=180,
                    plate=deep_plate,
                    col=col,
                    top_height=16,
                    bottom_height=2,
                    steps=2)
        p300_multi.drop_tip()

    protocol.delay(minutes=5, msg="Lysing cells...")

    # Transfer Neutralizing buffer
    protocol.comment("Dispensing neutralisation buffer.")

    for col in cols:
        buffer_well = 0
        p300_multi.pick_up_tip()
        transfer_large_vol(pipette=p300_multi,
                           total_vol=450,
                           source_well=neutralisation[buffer_well // 2],
                           dest_plate=deep_plate,
                           col=col)
        mix_up_down(pipette=p300_multi,
                    vol=180,
                    plate=deep_plate,
                    col=col,
                    top_height=18,
                    bottom_height=2,
                    steps=1)
        p300_multi.drop_tip()
        buffer_well += 1

    p300_multi.flow_rate.aspirate = 150

    # Clearing Beads
    protocol.comment("Mixing with Clearing Beads.")

    p300_multi.pick_up_tip()
    p300_multi.mix(5, 100, clearing_beads)  # Mix before dispensing beads

    for col in cols:
        check_tip_attached(p300_multi)
        p300_multi.mix(1, 100, clearing_beads)
        p300_multi.transfer(50, clearing_beads, deep_plate[col][0], new_tip='never')
        mix_up_down(pipette=p300_multi,
                    vol=180,
                    plate=deep_plate,
                    col=col,
                    top_height=18,
                    bottom_height=2,
                    steps=2)
        p300_multi.drop_tip()

    magmol.engage(height=mag_engage_height)
    protocol.delay(minutes=5, msg="Settling Clearing Beads...")

    # Transfer supernatant to second half of plate
    protocol.comment("Transferring supernatant.")

    p300_multi.well_bottom_clearance.dispense = 35

    for col in cols:
        p300_multi.pick_up_tip()
        p300_multi.default_speed = 200  # Set head speed to half of default to avoid dripping
        transfer_large_vol(pipette=p300_multi,
                           total_vol=750,
                           source_well=deep_plate[col][0].bottom(z=10),
                           dest_plate=deep_plate,
                           col=col+6,
                           speed_asp=30)
        # p300_multi.transfer([190, 190, 190, 180], deep_plate[col][0].bottom(z=10),
        #                     [deep_plate[c] for c in [col+6, col+6, col+6, col+6]], new_tip='never', air_gap=10)
        p300_multi.default_speed = 400  # Return head speed to default
        # p300_multi.blow_out(deep_plate[col+6][0])
        # p300_multi.touch_tip()
        p300_multi.drop_tip()

    magmol.disengage()

    # Redefine start and end columns in second half of plate
    start_col = start_col + 6
    end_col = end_col + 6
    cols = range(start_col, end_col)

    # Binding Beads
    protocol.comment("Mixing with Binding Beads.")

    p300_multi.pick_up_tip()
    p300_multi.mix(5, 100, binding_beads)  # Mix before dispensing beads

    for i in range(2):
        for col in cols:
            check_tip_attached(p300_multi)
            if i == 0:
                p300_multi.mix(1, 200, binding_beads)
                p300_multi.transfer(30, binding_beads, deep_plate[col][0], new_tip='never')
            mix_up_down(pipette=p300_multi,
                        vol=180,
                        plate=deep_plate,
                        col=col,
                        top_height=16,
                        bottom_height=2,
                        steps=2)
            p300_multi.drop_tip()

        protocol.delay(minutes=4, msg="Incubating with Binding Beads...")

    magmol.engage(height=mag_engage_height)
    protocol.delay(minutes=1, msg="Settling Binding Beads...")

    # Waste
    protocol.comment("Discarding supernatant.")

    for col in cols:
        discard_supernatant(p300_multi,
                           800,
                           deep_plate,
                           col,
                           waste,
                           bottom_offset=1)

    magmol.disengage()

    # Endo-Wash
    protocol.comment("Washing with Endo-Wash buffer.")
    for col in cols:
        p300_multi.pick_up_tip()
        p300_multi.transfer(200, endo_wash, deep_plate[col][0], new_tip='never')
        p300_multi.mix(5, 150, deep_plate[col][0].bottom(z=2), rate=2)
        p300_multi.drop_tip()

    magmol.engage(height=mag_engage_height)
    protocol.delay(minutes=1, msg="Beads settle...")

    # Waste
    protocol.comment("Discarding supernatant.")

    for col in cols:
        discard_supernatant(p300_multi,
                           220,
                           deep_plate,
                           col,
                           waste,
                           bottom_offset=1)

    magmol.disengage()

    # Wash 2x with Zyppy Wash
    wash_step = ["1st", "2nd"]

    for _ in range(2):
        protocol.comment("Washing %s time with Zyppy-Wash buffer." % wash_step[_])
        for col in cols:
            buffer_well = 0
            p300_multi.pick_up_tip()
            transfer_large_vol(pipette=p300_multi,
                               total_vol=380,
                               source_well=zyppy_wash[buffer_well],
                               dest_plate=deep_plate,
                               col=col)
            p300_multi.mix(5, 180, deep_plate[col][0].bottom(z=2), rate=2)
            p300_multi.drop_tip()
            buffer_well += 1

        magmol.engage(height=mag_engage_height)
        protocol.delay(minutes=1, msg="Beads settle...")

        # Waste
        protocol.comment("Discarding supernatant.")
        for col in cols:
            discard_supernatant(p300_multi,
                               400,
                               deep_plate,
                               col,
                               waste,
                               bottom_offset=1)

        magmol.disengage()

    protocol.delay(minutes=drying_time, msg="Drying beads for %d min. Protocol will resume automatically." % drying_time)

    # Elution
    protocol.comment("Adding elution buffer.")

    p300_multi.flow_rate.aspirate = 150
    p300_multi.well_bottom_clearance.dispense = 2

    for i in range(2):
        for col in cols:
            p300_multi.pick_up_tip()
            if i == 0:
                p300_multi.transfer(elution_buffer_volume + 5, elution_buffer, deep_plate[col][0], new_tip='never')
            p300_multi.mix(20, elution_buffer_volume - 15, deep_plate[col][0].bottom(z=1))
            p300_multi.blow_out(deep_plate[col][0].bottom(z=5))
            p300_multi.drop_tip()
        protocol.delay(minutes=5, msg="Incubating for 5 min with elution buffer")

    magmol.engage(height=mag_engage_height)
    protocol.delay(minutes=2, msg="Beads settle...")

    protocol.comment("Elution of bound DNA.")

    p300_multi.well_bottom_clearance.dispense = 10  # set dispense height to 1 cm above bottom
    p300_multi.flow_rate.aspirate = 30  # reduce flow-rate to 1/5 of default
    p300_multi.flow_rate.dispense = 50

    for col in cols:
        deep_plate_offset = odd_or_even(deep_plate, col)  # determine offset
        p300_multi.pick_up_tip()
        p300_multi.transfer(elution_buffer_volume, deep_plate_offset, elution_plate[elution_col][0],
                            new_tip='never')
        p300_multi.blow_out(elution_plate[elution_col][0].top(z=-3))
        p300_multi.drop_tip()

    # Done
    magmol.disengage()
    protocol.home()
    protocol.comment("Plasmid purification finished. Please take out the elution plate (slot 3).")
