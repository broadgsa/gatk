package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

public enum EventType {
    BASE_SUBSTITUTION(0, "M", "Base Substitution"),
    BASE_INSERTION(1, "I", "Base Insertion"),
    BASE_DELETION(2, "D", "Base Deletion");

    public final int index;
    private final String representation;
    private final String longRepresentation;

    private EventType(int index, String representation, String longRepresentation) {
        this.index = index;
        this.representation = representation;
        this.longRepresentation = longRepresentation;
    }

    public static EventType eventFrom(int index) {
        switch (index) {
            case 0:
                return BASE_SUBSTITUTION;
            case 1:
                return BASE_INSERTION;
            case 2:
                return BASE_DELETION;
            default:
                throw new ReviewedStingException(String.format("Event %d does not exist.", index));
        }        
    }
    
    public static EventType eventFrom(String event) {
        for (EventType eventType : EventType.values())
            if (eventType.representation.equals(event))
                return eventType;

        throw new ReviewedStingException(String.format("Event %s does not exist.", event));
    }

    @Override
    public String toString() {
        return representation;
    }

    public String prettyPrint() {
        return longRepresentation;
    }
}